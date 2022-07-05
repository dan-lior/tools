
#include <iostream>
#include <vector>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"

int main()
{
    std::cerr << "hello" << std::endl;

    // const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2_nondicom_data/");
    // std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> gwa_timeslices = read_nondicom(nondicom_directory);
    // GridWithAttribute<GridIso<3>, Vector<3>> grid = gwa_timeslices[8];

    const std::string filename("/home/dan/git-repos/tools/MRI/vtk_data/test.vtk");
    

    Index<3> n(2,2,2);
    std::vector<double> attribute({88,88,88,88,88,88,88,88});

    GridIso<3> grid(n);
    GridWithAttribute<GridIso<3>, double> gwa(grid);
    gwa.set_attribute(attribute);

    GridExport::export_to_vtk(gwa, filename);

    assert(Miscellaneous::file_exists(filename));

    std::cerr << "goodbye" << std::endl;
    return 0;
}