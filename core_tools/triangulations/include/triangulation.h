    #pragma once

    #include "common_defs.h"
    #include "disjoint_union.h"
    #include "misc.h"
    #include "affinity.h"

    class Polyline;

    struct Triangulation
    {
        // assumes that faces are all triangles
        void read_from_obj(const std::string& filename);

        void write_to_obj(const std::string& filename);

        Triangulation();
        Triangulation(const std::vector<Vector<3>>& vertices, const std::vector<Index<3>>& triangles, const std::vector<Vector<3>>& normals);

        Triangulation transform(const Affinity<3,3>& affinity) const;

        Triangulation extract_subtriangulation(const std::vector<uint64_t>& vertex_indices, uint64_t mode) const;

        // Vector<3> normal is a unit vector normal to the two slab faces
        Triangulation extract_slab(const Polyline& polyline, uint64_t iv, const Vector<3>& normal, uint64_t slab_radius, uint64_t mode) const;

        uint64_t closest_vertex(Vector<3> q) const;

        DisjointUnion connectivity() const;

        Triangulation extract_connected_component(const Vector<3>& q) const;

        // computes the average minimum edge length of a triangle
        double average_min_edge_length() const;

        std::vector<Vector<3>> vertices;
        std::vector<Vector<3>> normals;
        std::vector<Index<3>> triangles;
    };


