#include "voxel_box.h"

namespace
{
    // label is just a short descriptive string that can be displayed when visualizing with certain third party software (such as VisIt)
    template<typename T>
    void write_velocity_slice_to_vtk(
        const Index<3>& n, 
        const std::vector<Vector<3>>& points, 
        const std::vector<Vector<3>>& velocities, 
        const std::string& label, 
        const std::string& filename)
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
        file << "DIMENSIONS" << " " << n(0) << " " << n(1) << " " << n(2) << std::endl;
        file << "POINTS" << " " << n(0)*n(1)*n(2) << " " << "float" << std::endl;

        for (auto p : points)
            file << p[0] << " " << p[1] << " " << p[2] << std::endl;

        file << "POINT_DATA" << " " << n(0)*n(1)*n(2) << std::endl;
        file << "VECTORS" << " " << label << " " << vtk_data_type << std::endl;

        for (auto v : velocities)
            file << v[0] << " " << v[1] << " " << v[2] << std::endl;

        std::cerr << std::endl;
        file.close();

    }
}

template<typename GridType, typename T>
void GridExport::export_to_vtk(const GridWithAttribute<GridType, T>& gwa, const std::string& filename)
{
    std::cerr << "The function template export_to_vtk is currently only defined when its template parameters are: GridWithAttribute<GridIso<3>, Vector<3> " << std::endl;
}

template <>
void GridExport::export_to_vtk(const GridWithAttribute<GridIso<3>, Vector<3>>& gwa, const std::string& filename)
{
    const std::string& label("fluid velocity from 4dmri"); 

    std::vector<Vector<3>> points;
    for (uint64_t rank = 0; rank < gwa.grid.indexer.num_cells(); rank++)
    {
        Index<3> index = gwa.grid.indexer.to_index(rank);
        Vector<3> fractional_index = index_to_vector<3>(index);
        Vector<3> point = gwa.grid.position(fractional_index);
        points.push_back(point);
    }
 
    write_velocity_slice_to_vtk<double>(gwa.grid.indexer.n, points, gwa.get_attribute(), label, filename);
}
