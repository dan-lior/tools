#pragma once

#include "common_defs.h"
#include "nlohmann/json.hpp"

namespace Miscellaneous
{
    nlohmann::json read_metadata(const std::string& filename);

    // point clouds are always stored in vtk files as 3d points. 
    // here, 1d and 2d points are represented as 3d points on the x-axis and xy-plane, respectively

    template<uint64_t dim>
    std::vector<Vector<dim>> read_points_from_vtk_file(const std::string& filename_vtk);

    template<uint64_t dim>
    void write_points_to_vtk_file(const std::vector<Vector<dim>>& vertices, const std::string& filename_vtk, bool connect_the_dots = false);

    bool validate_filename_extension(const std::string& filename, const std::string& expected_extension);

    template<typename T>
    void read_file(std::vector<T>& store, const std::string& filename, uint64_t expected_number_of_elements);

    void write_to_file(const std::vector<std::vector<double>>& xss, const std::vector<std::string>& header, const std::string& filename);

    // copied from "map_masked_velocity_to_vtk.cpp", 3/18/2022
    template<typename T>
    void write_scalar_slice_to_vtk(uint64_t nx, uint64_t ny, uint64_t nz, const std::vector<std::array<double, 3>>& points, 
    const std::vector<T>& scalar_data, const std::string& label, const std::string& filename);

    bool file_exists(const std::string& filename);

    bool directory_exists(const std::string& directory_name);

    void verify_directory(const std::string& directory);
    void initialize_directory(const std::string& directory);
    void verify_file(const std::string& filename);

}


#include "misc.hpp"