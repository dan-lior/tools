#include "voxel_box.h"


template<typename GridType, typename T>
void GridExport::export_to_vtk(const GridWithAttribute<GridType, T>& grid, const std::string& filename)
{
    std::cerr << "defualt " << std::endl;
}

template <>
void GridExport::export_to_vtk(const GridWithAttribute<GridIso<3>, Vector<3>>& grid, const std::string& filename)
{
    std::cerr << "specialization" << std::endl;
}
