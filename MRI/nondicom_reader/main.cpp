
#include <iostream>
#include <vector>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"

int main()
{
    std::cerr << "hello" << std::endl;

    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2_nondicom_data/");
    std::vector<LabelledGrid<3, 3, Vector<3>>> gwa_timeslices = read_nondicom(nondicom_directory);

    const std::string vtk_directory("/home/dan/git-repos/tools/MRI/vtk_data/patient2/");

    for (uint64_t i=0; i<gwa_timeslices.size(); ++i)
    {
        std::string filename(vtk_directory + std::to_string(i) + std::string(".vtk"));
        std::cerr << "about to export timeslice " << i << std::endl;
        gwa_timeslices[i].export_to_vtk(filename);
        std::cerr << "finished exporting timeslice " << i  << std::endl;
    }

    std::cerr << "goodbye" << std::endl;
    return 0;
}