
#include <iostream>
#include <vector>
#include <cassert>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"
#include "math.h"
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
    void write_vector_slice_to_vtk(
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

    // during export, the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the attribute itself, however, is NOT transformed. 
    template<uint64_t dim_target, uint64_t dim_source, typename T>
    void export_labelled_grid_to_vtk(LabelledGrid<dim_target, dim_source, T> labelled_grid, const std::string& filename)
    {
//        static_assert(false, "shouldn't get here");
       static_assert(!(std::is_same<T, Vector<3>>::value || std::is_same<T, double>::value, "vtk export error: can only export double or Vector<3> attributes"));
//        static_assert(!(dim == 3 || dim == 2, "vtk export error: can only export dimension 2 or 3 objects"));
    }

    template<>
    void export_labelled_grid_to_vtk(LabelledGrid<3, 2, double> labelled_grid, const std::string& filename)
    {
        const Index<2> n = labelled_grid.indexer.n;
        Index<3> n2;
        n2(0) = n(0);
        n2(1) = n(1);
        n2(2) = 1;
        
        const Affinity<3,2> affinity = labelled_grid.affinity;
        const MatrixRect<3,2> A = affinity.linear_part;
        const Vector<3> b = affinity.translation_part;
        Matrix<3> A2 = Matrix<3>::Zero();
        A2(0,0) = A(0,0);
        A2(0,1) = A(0,1);
        A2(1,0) = A(1,0);
        A2(1,1) = A(1,1);
        A2(2,0) = A(2,0);
        A2(2,1) = A(2,1);

        Grid<3, 3> thickened_grid(n2, Affinity<3,3>(A2,b));
        LabelledGrid<3,3, double> labelled_thickened_grid(thickened_grid, labelled_grid.get_labels());

        std::cerr << "exporting scalar data" << std::endl;
        write_scalar_slice_to_vtk(labelled_thickened_grid.indexer.n, labelled_thickened_grid.all_positions(), labelled_thickened_grid.get_labels(), "some_scalar_quantity", filename);
    }


    template<>
    void export_labelled_grid_to_vtk(LabelledGrid<3, 2, Vector<3>> labelled_grid, const std::string& filename)
    {
        const Index<2> n = labelled_grid.indexer.n;
        Index<3> n2;
        n2(0) = n(0);
        n2(1) = n(1);
        n2(2) = 1;
        
        const Affinity<3,2> affinity = labelled_grid.affinity;
        const MatrixRect<3,2> A = affinity.linear_part;
        const Vector<3> b = affinity.translation_part;
        Matrix<3> A2 = Matrix<3>::Zero();
        A2(0,0) = A(0,0);
        A2(0,1) = A(0,1);
        A2(1,0) = A(1,0);
        A2(1,1) = A(1,1);
        A2(2,0) = A(2,0);
        A2(2,1) = A(2,1);

        Grid<3, 3> thickened_grid(n2, Affinity<3,3>(A2,b));
        LabelledGrid<3,3, Vector<3>> labelled_thickened_grid(thickened_grid, labelled_grid.get_labels());

        std::cerr << "exporting scalar data" << std::endl;
        write_vector_slice_to_vtk(labelled_thickened_grid.indexer.n, labelled_thickened_grid.all_positions(), labelled_thickened_grid.get_labels(), "some_vector_quantity", filename);
    }



    template<>
    void export_labelled_grid_to_vtk(LabelledGrid<3, 3, double> labelled_grid, const std::string& filename)
    {
        std::cerr << "exporting scalar data" << std::endl;
        write_scalar_slice_to_vtk(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_scalar_quantity", filename);
    }

    template<>
    void export_labelled_grid_to_vtk(LabelledGrid<3, 3, Vector<3>> labelled_grid, const std::string& filename)
    {
        std::cerr << "exporting vector data" << std::endl;
        write_vector_slice_to_vtk(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_vector_quantity", filename);
    }
}


int main()
{
    // converts nondicom timeslices to vtk timeslices

    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2/");
    const std::string vtk_directory("/home/dan/git-repos/tools/MRI/vtk_data/patient2/");
    const std::string vtk_slice_directory = vtk_directory + std::string("slices/");

    assert(Miscellaneous::directory_exists(nondicom_directory));
    assert(Miscellaneous::directory_exists(vtk_directory));
    assert(Miscellaneous::directory_exists(vtk_slice_directory));

    std::vector<LabelledGrid<3, 3, Vector<3>>> velocity_timeslices = read_nondicom(nondicom_directory);

    // for (uint64_t i=0; i<velocity_timeslices.size(); ++i)
    // {
    //     std::string filename(vtk_directory + std::to_string(i) + std::string(".vtk"));
    //     std::cerr << "about to export timeslice " << i << std::endl;
    //     velocity_timeslices[i].export_to_vtk(filename);
    //     std::cerr << "finished exporting timeslice " << i  << std::endl;
    // }

    LabelledGrid<3,3, Vector<3>> labelled_cube = velocity_timeslices[0];


    // MatrixRect<3,3> A;
    // A = MatrixRect<3,3>::Identity();
    // A(0,0) = 0.866;
    // A(0,1) = 0.002;
    // A(0,2) = 1.98;
    // A(1,0) = 1.8;
    // A(1,1) = -0.001;
    // A(1,2) = -0.95;
    // A(2,0) = 0;
    // A(2,1) = -2;
    // A(2,2) = 0.002;
    // Vector<3> b;
    // b(0) = -78.5;
    // b(1) = -48.5;
    // b(2) = 187.0;
    const Index<3> n = labelled_cube.indexer.n;
    const Affinity<3,3> dicom_transformation = labelled_cube.affinity;
    //const Affinity<3,3> dicom_transformation = Affinity<3,3>(A,b);
    // Grid<3,3> cube(n, dicom_transformation);



    typedef Vector<3> LabelType;
    // const LabelType default_value = 888.8;

    // LabelledGrid<3,3, LabelType> labelled_cube(cube, std::vector<LabelType>(cube.indexer.num_cells(), default_value));

    const Vector<3> position_of_origin = labelled_cube.position(index_to_vector<3>(Index<3>(0,0,0)));
    const Vector<3> position_of_antiorigin = labelled_cube.position(index_to_vector<3>(n - Index<3>(1,1,1)));
    const Vector<3> position_of_center = 0.5*(position_of_origin + position_of_antiorigin);

    std::cerr << "position_of_origin: \t" <<  position_of_origin.transpose() << std::endl;
    std::cerr << "position_of_antiorigin: \t" <<  position_of_antiorigin.transpose() << std::endl;
    std::cerr << "position_of_center: \t" <<  position_of_center.transpose() << std::endl;

    const Affinity<3,3> rotation_through_origin(MiscellaneousMath::rotation(Vector<3>(1,0,0), M_PI/6));
    const Affinity<3,3> translation(Matrix<3>::Identity(), position_of_center);
//    const Affinity<3,3> rotation = (translation.operator*<3>(rotation_through_origin)).operator*<3>(translation.inverse());
    const Affinity<3,3> rotation = rotation_through_origin;

    const Affinity<3,2> slice_transformation = (rotation.operator*<3>(dicom_transformation)).operator*<2>(Affinity<3,2>());
    
    const Grid<3,2> slice(Index<2>(n(0), n(1)), slice_transformation);
    

    auto rank_map = slice.rank_map(labelled_cube);

    assert(!rank_map.empty());

    std::vector<LabelType> labels(slice.indexer.num_cells(), LabelType());
    for (auto iter = rank_map.begin(); iter != rank_map.end(); ++iter)
        labels[iter->first] = -1*labelled_cube.labels[iter->second];

    LabelledGrid<3,2, LabelType> labelled_slice(slice, labels);


    export_labelled_grid_to_vtk(labelled_cube, vtk_slice_directory + std::string("cube.vtk"));
    export_labelled_grid_to_vtk(labelled_slice, vtk_slice_directory + std::string("slice.vtk"));


    return 0;
}