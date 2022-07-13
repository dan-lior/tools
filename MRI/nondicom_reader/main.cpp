
#include <iostream>
#include <vector>
#include <cassert>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"
//#include "math.h"

namespace
{

    Grid<3,2> generate_slice_from_cube_demo(const Grid<3,3>& cube)
    {
        const Index<3> n = cube.indexer.n;
        const Affinity<3,3> dicom_transformation = cube.affinity;

        const Vector<3> position_of_origin = cube.position(index_to_vector<3>(Index<3>(0,0,0)));
        const Vector<3> position_of_antiorigin = cube.position(index_to_vector<3>(n - Index<3>(1,1,1)));
        const Vector<3> position_of_center = 0.5*(position_of_origin + position_of_antiorigin);


        const Affinity<3,3> rotation_through_origin(MiscellaneousMath::rotation(Vector<3>(1,0,0), M_PI/6));
        const Affinity<3,3> translation(Matrix<3>::Identity(), position_of_center);
    //    const Affinity<3,3> rotation = (translation.operator*<3>(rotation_through_origin)).operator*<3>(translation.inverse());
        const Affinity<3,3> rotation = rotation_through_origin;

        const Affinity<3,2> slice_transformation = (rotation.operator*<3>(dicom_transformation)).operator*<2>(Affinity<3,2>());
        
        const Grid<3,2> slice(Index<2>(n(0), n(1)), slice_transformation);

        return slice;
    }

}

int main()
{
    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2/");
    const std::string vtk_directory("/home/dan/git-repos/tools/MRI/vtk_data/patient2/");
    const std::string vtk_slice_directory = vtk_directory + std::string("slices/");

    assert(Miscellaneous::directory_exists(nondicom_directory));
    assert(Miscellaneous::directory_exists(vtk_directory));
    assert(Miscellaneous::directory_exists(vtk_slice_directory));

    std::vector<LabelledGrid<3, 3, Vector<3>>> velocity_timeslices = read_nondicom(nondicom_directory);

    // choose a single timeslice from which to determine geometry for the space slice

    Grid<3,2> slice = generate_slice_from_cube_demo(velocity_timeslices[0]);

    auto rank_map = slice.rank_map(velocity_timeslices[0]);

    assert(!rank_map.empty());

    // use the rank map to generate space slice slices for each time slice
    // export both the cube and the slice to vtk file

    for (uint64_t i=0; i<velocity_timeslices.size(); ++i)
    {
        LabelledGrid<3,3, Vector<3>> labelled_cube = velocity_timeslices[i];

        auto cube_labels = labelled_cube.get_labels();

        std::vector<Vector<3>> labels(slice.indexer.num_cells(), Vector<3>());
        for (auto iter = rank_map.begin(); iter != rank_map.end(); ++iter)
            labels[iter->first] = cube_labels[iter->second];
        LabelledGrid<3,2, Vector<3>> labelled_slice(slice, labels);

        labelled_cube.export_to_vtk(vtk_directory + std::to_string(i) + std::string(".vtk"));
        labelled_slice.export_to_vtk(vtk_slice_directory + std::to_string(i) + std::string(".vtk"));

        std::cerr << "finished timeslice " << i << " of " << velocity_timeslices.size() << std::endl;
    }

    return 0;
}