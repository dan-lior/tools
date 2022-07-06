
#include <iostream>
#include <vector>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"

int main()
{
    // converts nondicom timeslices to vtk timeslices

    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2/");
    const std::string vtk_directory("/home/dan/git-repos/tools/MRI/vtk_data/patient2/");

    std::vector<LabelledGrid<3, 3, Vector<3>>> velocity_timeslices = read_nondicom(nondicom_directory);

    for (uint64_t i=0; i<velocity_timeslices.size(); ++i)
    {
        std::string filename(vtk_directory + std::to_string(i) + std::string(".vtk"));
        std::cerr << "about to export timeslice " << i << std::endl;
        velocity_timeslices[i].export_to_vtk(filename);
        std::cerr << "finished exporting timeslice " << i  << std::endl;
    }

    return 0;
}