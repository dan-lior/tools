
#include <iostream>
#include <vector>
#include "voxel_box.h"
#include "read_nondicom.h"
#include "misc.h"

int main()
{
    std::cerr << "hello" << std::endl;

    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient2_nondicom_data/");
    std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> gwa_timeslices = read_nondicom(nondicom_directory);

    const std::string vtk_directory("/home/dan/git-repos/tools/MRI/vtk_data/");

    for (uint64_t i=0; i<gwa_timeslices.size(); ++i)
    {
        std::string filename(vtk_directory + std::to_string(i) + std::string(".vtk"));
        GridExport::export_to_vtk(gwa_timeslices[i], filename);
    }

    // Index<3> n(2,2,2);
    // GridIso<3> grid(n);

    // GridWithAttribute<GridIso<3>, double> gwsa(grid);
    // GridWithAttribute<GridIso<3>, Vector<3>> gwva(grid);

    // std::vector<double> scalar_attribute({89,90,91,92,93,94,95,96});
    // std::vector<Vector<3>> vector_attribute({
    //     Vector<3>(2,3,4), Vector<3>(2,3,4), Vector<3>(2,3,4), Vector<3>(2,3,4), 
    //     Vector<3>(2,3,4), Vector<3>(2,0,4), Vector<3>(0,3,4), Vector<3>(2,3,0) 
    //     });
    
    // gwsa.set_attribute(scalar_attribute);
    // gwva.set_attribute(vector_attribute);

    // GridExport::export_to_vtk(gwsa, std::string("/home/dan/git-repos/tools/MRI/vtk_data/test_scalar.vtk"));
    // GridExport::export_to_vtk(gwva, std::string("/home/dan/git-repos/tools/MRI/vtk_data/test_vector.vtk"));

    std::cerr << "goodbye" << std::endl;
    return 0;
}