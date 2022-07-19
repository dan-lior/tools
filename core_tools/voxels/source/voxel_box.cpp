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

template<>
void LabelledGrid<3,3,double>::export_labelled_grid_to_vtk(const LabelledGrid<3,3,double>& labelled_grid, const std::string& filename)
{
    std::cerr << "exporting scalar data" << std::endl;

    Miscellaneous::write_scalar_slice_to_vtk<double>(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_scalar_quantity", filename);
}

template<>
void LabelledGrid<3,3,uint64_t>::export_labelled_grid_to_vtk(const LabelledGrid<3,3,uint64_t>& labelled_grid, const std::string& filename)
{
    std::cerr << "exporting scalar data" << std::endl;

    Miscellaneous::write_scalar_slice_to_vtk<uint64_t>(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_scalar_quantity", filename);
}

template<>
void LabelledGrid<3,3,Vector<3>>::export_labelled_grid_to_vtk(const LabelledGrid<3,3,Vector<3>>& labelled_grid, const std::string& filename)
{
    std::cerr << "exporting vector data" << std::endl;
    Miscellaneous::write_vector_slice_to_vtk(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_vector_quantity", filename);
}


// template<>
// void LabelledGrid<3, 2, double>::export_labelled_grid_to_vtk(const std::string& filename) const
// {
//     const Index<2> n = indexer.n;
//     Index<3> n2;
//     n2(0) = n(0);
//     n2(1) = n(1);
//     n2(2) = 1;
    
//     const Affinity<3,2> affinity = affinity;
//     const MatrixRect<3,2> A = affinity.linear_part;
//     const Vector<3> b = affinity.translation_part;
//     Matrix<3> A2 = Matrix<3>::Zero();
//     A2(0,0) = A(0,0);
//     A2(0,1) = A(0,1);
//     A2(1,0) = A(1,0);
//     A2(1,1) = A(1,1);
//     A2(2,0) = A(2,0);
//     A2(2,1) = A(2,1);

//     Grid<3, 3> thickened_grid(n2, Affinity<3,3>(A2,b));
//     LabelledGrid<3,3, double> labelled_thickened_grid(thickened_grid, get_labels());

//     std::cerr << "exporting scalar data" << std::endl;
//     Miscellaneous::write_scalar_slice_to_vtk(labelled_thickened_grid.indexer.n, labelled_thickened_grid.all_positions(), labelled_thickened_grid.get_labels(), "some_scalar_quantity", filename);
// }

// template<>
// void LabelledGrid<3, 2, Vector<3>>::export_labelled_grid_to_vtk(const std::string& filename) const
// {
//     const Index<2> n = indexer.n;
//     Index<3> n2;
//     n2(0) = n(0);
//     n2(1) = n(1);
//     n2(2) = 1;
    
//     const Affinity<3,2> affinity = affinity;
//     const MatrixRect<3,2> A = affinity.linear_part;
//     const Vector<3> b = affinity.translation_part;
//     Matrix<3> A2 = Matrix<3>::Zero();
//     A2(0,0) = A(0,0);
//     A2(0,1) = A(0,1);
//     A2(1,0) = A(1,0);
//     A2(1,1) = A(1,1);
//     A2(2,0) = A(2,0);
//     A2(2,1) = A(2,1);

//     Grid<3, 3> thickened_grid(n2, Affinity<3,3>(A2,b));
//     LabelledGrid<3,3, Vector<3>> labelled_thickened_grid(thickened_grid, get_labels());

//     std::cerr << "exporting vector data" << std::endl;
//     Miscellaneous::write_vector_slice_to_vtk(labelled_thickened_grid.indexer.n, labelled_thickened_grid.all_positions(), labelled_thickened_grid.get_labels(), "some_vector_quantity", filename);
// }



