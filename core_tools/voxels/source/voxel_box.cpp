#include "voxel_box.h"

namespace
{
    void write_common_vtk_part(std::ofstream& file, const Index<3>& n, const std::vector<Vector<3>>& points)
    {
        file << "# vtk DataFile Version 2.3" << std::endl; // prior to this version, numComp was not supported
        file << "DANWUZHERE" << std::endl;
        file << "ASCII" << std::endl;
        file << "DATASET STRUCTURED_GRID" << std::endl;
        file << "DIMENSIONS" << " " << n(0) << " " << n(1) << " " << n(2) << std::endl;
        file << "POINTS" << " " << n(0)*n(1)*n(2) << " " << "double" << std::endl;

        for (auto p : points)
            file << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }

    // label is just a short descriptive string that is displayed when visualizing vtk files with certain third party software (e.g. VisIt)
    void write_scalar_slice_to_vtk(
        const Index<3>& n, 
        const std::vector<Vector<3>>& points, 
        const std::vector<double>& scalar_data, 
        const std::string& label, 
        const std::string& filename)
    {
        std::ofstream file(filename, std::ios::binary);
        write_common_vtk_part(file, n, points);        

        file << "POINT_DATA " << n(0)*n(1)*n(2) << std::endl;
        file << "SCALARS " << label << " double" << " " << 1 << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        for (auto v : scalar_data) 
            file << v << std::endl;

        file.close();
    }

    // label is just a short descriptive string that is displayed when visualizing vtk files with certain third party software (e.g. VisIt)
    void write_velocity_slice_to_vtk(
        const Index<3>& n, 
        const std::vector<Vector<3>>& points, 
        const std::vector<Vector<3>>& velocities, 
        const std::string& label, 
        const std::string& filename)
    {
        std::ofstream file(filename, std::ios::binary);
        write_common_vtk_part(file, n, points);        

        file << "POINT_DATA " << n(0)*n(1)*n(2) << std::endl;
        file << "VECTORS " << label << " double" << std::endl;

        for (auto v : velocities)
            file << v[0] << " " << v[1] << " " << v[2] << std::endl;

        file.close();
    }



}




template<typename GridType, typename T>
void GridExport::export_to_vtk(const GridWithAttribute<GridType, T>& gwa, const std::string& filename)
{
    std::cerr << "The function template export_to_vtk is currently only defined when its template parameters are either: " << std::endl;
    std::cerr << "GridWithAttribute<GridIso<3>, Vector<3> " << std::endl;
    std::cerr << "or " << std::endl;
    std::cerr << "GridWithAttribute<GridIso<3>, double " << std::endl;
}

template <>
void GridExport::export_to_vtk(const GridWithAttribute<GridIso<3>, double>& gwa, const std::string& filename)
{
    const std::string& label("some_scalar_quantity"); 
    write_scalar_slice_to_vtk(gwa.grid.indexer.n, gwa.grid.all_positions(), gwa.get_attribute(), label, filename);
}

template <>
void GridExport::export_to_vtk(const GridWithAttribute<GridIso<3>, Vector<3>>& gwa, const std::string& filename)
{
    const std::string& label("some_vector_quantity"); 
    write_velocity_slice_to_vtk(gwa.grid.indexer.n, gwa.grid.all_positions(), gwa.get_attribute(), label, filename);
}

