
// copied from "map_masked_velocity_to_vtk.cpp", 3/18/2022
template<typename T>
void Miscellaneous::write_scalar_slice_to_vtk(uint64_t nx, uint64_t ny, uint64_t nz, const std::vector<std::array<double, 3>>& points, 
const std::vector<T>& scalar_data, const std::string& label, const std::string& filename)
{
    std::string vtk_data_type;
    if (typeid(T) == typeid(bool)) vtk_data_type = "bit";
    else if (typeid(T) == typeid(double)) vtk_data_type = "double";
    else if (typeid(T) == typeid(float)) vtk_data_type = "float";
    else if (typeid(T) == typeid(char)) vtk_data_type = "char";
    else if (typeid(T) == typeid(short)) vtk_data_type = "short";
    else if (typeid(T) == typeid(int)) vtk_data_type = "int";
    else if (typeid(T) == typeid(long)) vtk_data_type = "long";
    else if (typeid(T) == typeid(unsigned char)) vtk_data_type = "unsigned_char";
    else if (typeid(T) == typeid(unsigned short)) vtk_data_type = "unsigned_short";
    else if (typeid(T) == typeid(unsigned int)) vtk_data_type = "unsigned_int";
    else if (typeid(T) == typeid(unsigned long)) vtk_data_type = "unsigned_long";
    else
    {
        std::cerr << "error: unsupported type" << std::endl;	
        exit(1);	
    }

    std::ofstream file(filename, std::ios::binary);
    file << "# vtk DataFile Version 2.3" << std::endl; // prior to this version, numComp was not supported
    file << "DANWUZHERE" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_GRID" << std::endl;
    file << "DIMENSIONS" << " " << nx << " " << ny << " " << nz << std::endl;
    file << "POINTS" << " " << nx*ny*nz << " " << "double" << std::endl;

    for (auto p : points)
        file << p[0] << " " << p[1] << " " << p[2] << std::endl;

    file << "POINT_DATA" << " " << nx*ny*nz << std::endl;
    file << "SCALARS " << label << " " << vtk_data_type << " " << 1 << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    for (auto v : scalar_data) 
        file << v << std::endl;

    file.close();
}

template<uint64_t dim>
void Miscellaneous::write_points_to_vtk_file(const std::vector<Vector<dim>>& vertices, const std::string& filename_vtk, bool connect_the_dots)
{
    assert(Miscellaneous::validate_filename_extension(filename_vtk, "vtk"));

    std::ofstream file(filename_vtk);
    if (!file.is_open())
    {
        std::cerr << "file error with vtk file: " << filename_vtk << std::endl;
        exit(1);
    }

    file << "# vtk DataFile Version 4.2" << std::endl;
    file << "definition precedes existence" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET POLYDATA" << std::endl;
    file << "POINTS " << vertices.size() << " " << "float" << std::endl;

    switch (dim)
    {
    case 1:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  0 << " " <<  0 << " " <<  std::endl;
        break;
    
    case 2:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  vertex(1) << " " <<  0 << " " <<  std::endl;
        break;
    
    case 3:
        for (auto vertex : vertices)
            file << vertex(0) << " " <<  vertex(1) << " " <<  vertex(2) << " " <<  std::endl;
        break;
    
    default:
        std::cerr << "point dimensions other than 1,2,3 are currently not supported" <<std::endl;
        exit(1);
        break;
    }

    if (connect_the_dots)
    {
        file << "POLYGONS 1 " << vertices.size() << std::endl;
        for (uint64_t i=0; i<vertices.size(); ++i)
            file << vertices.size()-i << ' ';

        // "close the loop"
        file << vertices.size() << std::endl;
        file << std::endl;
    }

    file.close();
}

template<uint64_t dim>
std::vector<Vector<dim>> Miscellaneous::read_points_from_vtk_file(const std::string& filename_vtk)
{
    if (dim != 1 && dim != 2 && dim != 3)
    {
        std::cerr << "error: in file: " << filename_vtk << ", dim = " << dim << ", but point dimensions other than 1,2,3 are currently not supported" << std::endl;
        exit(1);
    }

    assert(Miscellaneous::validate_filename_extension(filename_vtk, "vtk"));

    std::ifstream file(filename_vtk);
    assert(file.is_open());

    // skip to the POINTS section
    std::string line;
    std::string first_word;
    std::string second_word;
    std::string third_word;
    do {
        std::getline(file, line);
        std::stringstream ss(line);
        ss >> first_word >> second_word >> third_word;
    } while (first_word.compare(std::string("POINTS")) != 0);
    uint64_t n = std::stoul(second_word);
    assert(third_word.compare("float") == 0);

    std::vector<Vector<dim>> vertices;

    if (dim == 3)
    {
        for (uint64_t i=0; i<n; ++i)
        {
            Vector<dim> v;
            for (uint64_t j=0; j<3; ++j)
            {
                std::string word;
                file >> word;
                float f = std::stof(word);
                v(j) = static_cast<double>(f);
            }
            vertices.push_back(v);
        }
    }
    else
    {
        for (uint64_t i=0; i<n; ++i)
        {
            Vector<dim> v;
            for (uint64_t j=0; j<3; ++j)
            {
                std::string word;
                file >> word;
                float f = std::stof(word);
                if (j<dim)
                {
                    if (f != 0.0)
                    {
                        std::cerr << "error with file: " << filename_vtk << std::endl;
                        std::cerr << "expected a file with data of dimension : " << dim << std::endl;
                        std::cerr << "component " << j << " was expected to be zero but read " << f << std::endl;
                        exit(1);
                    }
                    v(j) = static_cast<double>(f);
                }
            }
            vertices.push_back(v);
        }
    }

    file.close();

    return vertices;
}
