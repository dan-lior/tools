
#include "common_defs.h"
#include "voxel_box.h"
#include "nlohmann/json.hpp"

Index<3> grid_dimensions_from_nondicom_metadata(const nlohmann::json& nondicom_metadata)
{
    const uint64_t nx = nondicom_metadata["width"];
    const uint64_t ny = nondicom_metadata["height"];
    const uint64_t nz = nondicom_metadata["depth"];

    return Index<3>(nx,ny,nz);
}

AffinityIso<3> affinity_from_nondicom_metadata(const nlohmann::json& nondicom_metadata)
{
    // define a translation from the dicom data:
    Vector<3> b;
    b(0) = nondicom_metadata["image_position"]["x"];
    b(1) = nondicom_metadata["image_position"]["y"];
    b(2) = nondicom_metadata["image_position"]["z"];

    // define a rotation from the dicom data:
    // according to the dicom spec, the first two columns are orthogonal and of unit length
    Matrix<3> a;

    a(0,0) = nondicom_metadata["image_orientation"]["Xx"];
    a(1,0) = nondicom_metadata["image_orientation"]["Xy"];
    a(2,0) = nondicom_metadata["image_orientation"]["Xz"];

    a(0,1) = nondicom_metadata["image_orientation"]["Yx"];
    a(1,1) = nondicom_metadata["image_orientation"]["Yy"];
    a(2,1) = nondicom_metadata["image_orientation"]["Yz"];

    a.col(2) = a.col(0).cross(a.col(1));

    // define a scaling from the dicom data
    const double dx = nondicom_metadata["pixel_spacing"]["x"];
    const double dy = nondicom_metadata["pixel_spacing"]["y"];
    const double dz = nondicom_metadata["slice_thickness"];

    Matrix<3> d = Matrix<3>::Zero();
    d(0,0) = dx;
    d(1,1) = dy;
    d(2,2) = dz * (-1.0); // A hack to match conventions with 3d slicer

    return AffinityIso<3>(a*d, b);
}

GridIso<3> grid_from_nondicom_metadata(const nlohmann::json& nondicom_metadata)
{
    return GridIso<3>(grid_dimensions_from_nondicom_metadata(nondicom_metadata), affinity_from_nondicom_metadata(nondicom_metadata));
}

// GridWithAttribute<GridIso<3>, Vector<3>> read_nondicom(const std::string& nondicom_directory)
// {
//     //todo
// }


int main()
{

    std::cerr << "hello" << std::endl;
    return 0;
}