
#include <iostream>
#include <vector>
#include <cassert>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"

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

    // during export, the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the attribute itself, however, is NOT transformed. 
    template<uint64_t dim, typename T>
    void export_labelled_grid_to_vtk(LabelledGrid<dim, dim, T> labelled_grid, const std::string& filename)
    {
        static_assert(dim != 3 ||  (!std::is_same<T, Vector<3>>::value && !std::is_same<T, double>::value), "logic_error");
        std::cerr << "invalid template parameters for export" << std::endl;
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
        write_velocity_slice_to_vtk(labelled_grid.indexer.n, labelled_grid.all_positions(), labelled_grid.get_labels(), "some_vector_quantity", filename);
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

    // test slicing
    auto cube = velocity_timeslices[0];
    cube.labels = std::vector<Vector<3>>(cube.indexer.num_cells(), Vector<3>());
    
    const Vector<3> position_of_center_voxel = cube.position(index_to_vector<3>(cube.indexer.n / 2));

    const MatrixRect<3,2> xy_plane = MatrixRect<3,2>::Identity();
    const Affinity<3,2> affinity(xy_plane, Vector<3>(0,0,position_of_center_voxel(2)));
    const Index<2> n(200,200);
    const Grid<3,2> slice(n, affinity);

    auto rank_map = slice.rank_map(cube);

    for (auto iter = rank_map.begin(); iter != rank_map.end(); ++iter)
    {
        uint64_t other_rank = iter->second;
        cube.labels[other_rank] = 888*Vector<3>(1,0,0);
    }

//     for (uint64_t rank = 0; rank<cube.indexer.num_cells(); ++rank)
//     {
//         Index<3> index = cube.indexer.to_index(rank);
//         uint64_t ix = index(0);
//         uint64_t iy = index(1);
//         uint64_t iz = index(2);

//         Vector<3> fractional_index = index_to_vector<3>(index) + Vector<3>(0.5, 0.5, 0.5);
//         Vector<3> pos = cube.position(fractional_index);
//         double x = pos(0);
//         double y = pos(1);
//         double z = pos(2);

// //        if (nx/2 < ix && ix < nx/2+5)
//         if (-2 < y && y < 2)
//             cube.labels[rank] = 88*Vector<3>(1,0,0);
//         else
//             cube.labels[rank] = 0*Vector<3>(1,0,0);
//     }

	// const Affinity<3,2> sampling_plane(MatrixRect<3,2>(), Vector<3>(30,50,70));

	// const uint64_t nx = 100;
	// const uint64_t ny = 200;
	// const Index<2> n({nx,ny});
	// Grid<3,2> probe(n, sampling_plane);

    // LabelledGrid<3, 3, Vector<3>> slice = velocity_cube.clear_labels_outside_probe(probe);

    std::string slice_filename(vtk_slice_directory + std::string("slice.vtk"));

    // export_labelled_grid_to_vtk(slice, slice_filename);
    export_labelled_grid_to_vtk(cube, slice_filename);

    return 0;
}