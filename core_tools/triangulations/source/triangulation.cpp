#include "triangulation.h"
#include "polyline.h"
#include "disjoint_union.h"

void Triangulation::read_from_obj(const std::string& filename)
{
    assert(Miscellaneous::validate_filename_extension(filename, "obj"));
    std::ifstream file(filename);
    assert(file.is_open());

    std::string line;
    std::string first_word;

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string first_word;
        ss >> first_word;
        if (first_word.compare("v") == 0) // vertex coordinates 
        {
            Vector<3> vertex;
            for (unsigned int i=0; i<3; ++i)
            {
                std::string word;
                ss >> word;
                vertex(i) = std::stod(word);
            }            
            vertices.push_back(vertex);
        }
        else if (first_word.compare("vn") == 0) // vertex normal
        {
            Vector<3> normal;
            for (unsigned int i=0; i<3; ++i)
            {
                std::string word;
                ss >> word;
                normal(i) = std::stod(word);
            }            
            normals.push_back(normal);
        } 
        else if (first_word.compare("f") == 0) // ranks (starting at 1 of corner vertices 
        {
            Index<3> triangle;
            std::string word;
            for (unsigned int i=0; i<3; ++i)
            {
                ss >> word;
                size_t pos = word.find_first_of("/");
                triangle(i) = std::stoi(word.substr(0, pos)) - 1;
            }            
            triangles.push_back(triangle);
        }
    }
    file.close();

    // some validation
    assert(vertices.size() == normals.size());
    for (auto triangle : triangles)
    {
        assert(static_cast<uint64_t>(triangle.maxCoeff()) < vertices.size());
        uint64_t v0 = triangle(0);
        uint64_t v1 = triangle(1);
        uint64_t v2 = triangle(2);
        if (v0 == v1 || v1 == v2 || v1 == v0)
        {
            std::cerr << "improper triangle" << std::endl;
        }
        assert(v0 != v1);
        assert(v0 != v2);
        assert(v1 != v2);
    }
}

void Triangulation::write_to_obj(const std::string& filename)
{
    assert(Miscellaneous::validate_filename_extension(filename, "obj"));

    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "couldn't open file: " << filename << std::endl; 
        exit(1);
    }

    for (uint64_t i =0; i<vertices.size(); ++i)
    {
        file << "v " << vertices[i].transpose() << std::endl;
        file << "vn " << normals[i].transpose() << std::endl;
    }
    for (auto triangle : triangles) 
    {
        // .obj format uses 1-based indexing for faces
        file << "f " << 1+triangle(0) << " " << 1+triangle(1) << " " << 1+triangle(2) << std::endl;
    }

    file.close();
}

Triangulation::Triangulation() : 
    vertices(std::vector<Vector<3>>()), triangles(std::vector<Index<3>>()), normals(std::vector<Vector<3>>())
    {}

Triangulation::Triangulation(const std::vector<Vector<3>>& vertices_, const std::vector<Index<3>>& triangles_, const std::vector<Vector<3>>& normals_) : 
    vertices(vertices_), triangles(triangles_), normals(normals_)
    {}

Triangulation Triangulation::transform(const Affinity<3,3>& A) const
{
    std::vector<Vector<3>> vertices_transformed;
    for (auto v:vertices)
        vertices_transformed.push_back(A(v));

    std::vector<Index<3>> triangles_transformed = triangles;

    // apply only the linear part of the affinity to the normals of the triangulation
    std::vector<Vector<3>> normals_transformed;
    const Vector<3> translation = A(Vector<3>(0.0,0.0,0.0));
    for (auto v:normals)
        normals_transformed.push_back(A(v) - translation);

    return Triangulation(vertices_transformed, triangles_transformed, normals_transformed);
}


Triangulation Triangulation::extract_subtriangulation(const std::vector<uint64_t>& vertex_indices, uint64_t mode) const
{
    assert(!vertex_indices.empty());
    auto max_index = *std::max_element(vertex_indices.begin(), vertex_indices.end());
    if (max_index >= vertices.size())
    {
        std::cerr << "logic error" << std::endl;
    }

    std::set<uint64_t> indices_set(vertex_indices.begin(), vertex_indices.end());
    std::set<uint64_t> additional_indices;

    std::vector<Index<3>> selected_triangles;

    for (auto triangle : triangles)
    {
        uint64_t n=0;
        if (indices_set.find(triangle(0)) != indices_set.end()) n++;
        if (indices_set.find(triangle(1)) != indices_set.end()) n++;
        if (indices_set.find(triangle(2)) != indices_set.end()) n++;
        if (n >= mode) 
        {
            additional_indices.insert(triangle(0));
            additional_indices.insert(triangle(1));
            additional_indices.insert(triangle(2));                
            selected_triangles.push_back(triangle);
        }
    }

    indices_set.insert(additional_indices.begin(), additional_indices.end());
    
    std::map<uint64_t, uint64_t> dictionary;
    uint64_t rank=0;
    for (auto index : indices_set)
        dictionary[index] = rank++;
    assert(dictionary.size() == indices_set.size());


    std::vector<Vector<3>> vertices_rtn;
    std::vector<Vector<3>> normals_rtn;
    std::vector<Index<3>> triangles_rtn;

    for (auto index : indices_set)
    {
        vertices_rtn.push_back(vertices[index]);
        normals_rtn.push_back(normals[index]);
    }

    for (auto triangle : selected_triangles)
    {
        Index<3> reindexed_triangle;
        reindexed_triangle(0) = dictionary[triangle(0)];
        reindexed_triangle(1) = dictionary[triangle(1)];
        reindexed_triangle(2) = dictionary[triangle(2)];
        triangles_rtn.push_back(reindexed_triangle);
    }

    return Triangulation(vertices_rtn, triangles_rtn, normals_rtn);
}

Triangulation Triangulation::extract_slab(const Polyline& polyline, uint64_t iv, const Vector<3>& normal, uint64_t slab_radius, uint64_t mode) const
{
    const uint64_t i_first = (slab_radius > iv ? 0 : iv-slab_radius);
    const uint64_t i_last = std::min<uint64_t>(polyline.vertices.size()-1, iv + slab_radius);

    Vector<3> v_beg = polyline.vertices[i_first]; // a point on the first slab face
    Vector<3> v_end = polyline.vertices[i_last]; // a point on the second slab face

    std::vector<uint64_t> inside_indices;
    for (uint64_t i = 0; i<vertices.size(); ++i)
    {
        const Vector<3> p = vertices[i];
        if ((p-v_beg).dot(normal) * (p-v_end).dot(normal) < 0) // p lies on different sides of the two slab faces
            inside_indices.push_back(i);
    }

    return extract_subtriangulation(inside_indices, mode);
}

uint64_t Triangulation::closest_vertex(Vector<3> q) const
{
    double min_distance_sq = std::numeric_limits<double>::max();
    uint64_t i_min = 0;
    for(uint64_t i=0; i<vertices.size(); ++i)
    {
        double dist_sq = (q-vertices[i]).squaredNorm();
        if (dist_sq < min_distance_sq)
        {
            i_min = i;
            min_distance_sq = dist_sq;
        }
    }
    return i_min;
}

double Triangulation::average_min_edge_length() const
{
    if (triangles.empty())
    {
        std::cerr << "logic error" <<std::endl;
        exit(1);
    }

    double sum = 0;
    for (auto triangle : triangles)
    {
        Vector<3> v0 = vertices[triangle(0)];
        Vector<3> v1 = vertices[triangle(1)];
        Vector<3> v2 = vertices[triangle(2)];

        const double edge_length0 = (v1-v2).norm();
        const double edge_length1 = (v0-v2).norm();
        const double edge_length2 = (v1-v0).norm();

        sum += std::min<double>(edge_length0, std::min<double>(edge_length1, edge_length2));
    }

    return sum / triangles.size();
}



DisjointUnion Triangulation::connectivity() const
{
    DisjointUnion du(vertices.size());
    for (auto triangle : triangles)
    {
        uint64_t i = triangle(0);
        uint64_t j = triangle(1);
        uint64_t k = triangle(2);

        du.merge(i,j);
        du.merge(i,k);
    }
    du.cleanup();

    return du;
}

Triangulation Triangulation::extract_connected_component(const Vector<3>& q) const
{
    DisjointUnion connectivity1 = connectivity();
    uint64_t num_components = connectivity1.number_of_disjoint_subsets();
    assert(num_components != 0);

    Triangulation connected_component = *this;

    if (num_components > 1)
    {
        uint64_t i_min = closest_vertex(q);
        std::vector<uint64_t> component_vertices = connectivity1.disjoint_subset(i_min);
        connected_component = extract_subtriangulation(component_vertices, 3);

        // validation
        {
            DisjointUnion connectivity1 = connected_component.connectivity();
            uint64_t num_components = connectivity1.number_of_disjoint_subsets();
            if (num_components != 1)
                std::cerr << "logic error" << std::endl;
        }
    }

    assert(!connected_component.vertices.empty());

    return connected_component;
}
