
//#include <iostream>
//#include <vector>
//#include "voxel_box.h"
#include "read_nondicom.h"

int main()
{
    std::cerr << "hello" << std::endl;

    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2_nondicom_data/");
    std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> gwa_timeslices = read_nondicom(nondicom_directory);

    std::cerr << "goodbye" << std::endl;
    return 0;
}