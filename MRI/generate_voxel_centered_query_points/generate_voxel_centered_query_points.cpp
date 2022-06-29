// reads 4dmri data from file (as extracted from dicom files by my matlab code)
// reads a closed surface stl (segmented aorta)
// creates a voxel model, of the same dimensions as the 4dmri data, 
// populated with boolean values indicating which voxels are in the interior of the surface

#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <filesystem> 
#include "/home/dan/git-repos/json/single_include/nlohmann/json.hpp"
//#include "json.hpp"
using json = nlohmann::json;



json read_metadata(const std::string& filename)
{
	std::stringstream ss;
	std::ifstream in(filename);
	ss << in.rdbuf();
	in.close();
	return json::parse(ss.str());
}

const std::array<double, 3> apply_affine_transform(const std::array<std::array<double, 3> ,3>& a, const std::array<double, 3>& b, const std::array<double, 3>& x)
{
	std::array<double, 3> rtn;

	rtn[0] = b[0] + a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2];
	rtn[1] = b[1] + a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2];
	rtn[2] = b[2] + a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2];

	return rtn;
}

const std::array<double, 3> apply_dicom_transform(const json& metadata, const std::array<double, 3>& i, bool without_translation)
{
	const double dx = metadata["pixel_spacing"]["x"];
	const double dy = metadata["pixel_spacing"]["y"];
	const double dz = metadata["slice_thickness"];

	std::array<double, 3> x;
	x[0] = i[0]*dx;
	x[1] = i[1]*dy;
	x[2] = i[2]*dz * (-1.0); // A hack to match conventions with 3d slicer

	// define a translation from the dicom data:
	std::array<double, 3> b({0,0,0});
	if (!without_translation)
	{
		b[0] = metadata["image_position"]["x"];
		b[1] = metadata["image_position"]["y"];
		b[2] = metadata["image_position"]["z"];
	}

	// define a rotation from the dicom data:
	// according to the dicom convention, these first two columns are orthogonal and of unit length
	// to get an rotation matrix, 3rd column is the cross product of the first two
	std::array<std::array<double, 3> ,3> a;

	a[0][0] = metadata["image_orientation"]["Xx"];
	a[1][0] = metadata["image_orientation"]["Xy"];
	a[2][0] = metadata["image_orientation"]["Xz"];

	a[0][1] = metadata["image_orientation"]["Yx"];
	a[1][1] = metadata["image_orientation"]["Yy"];
	a[2][1] = metadata["image_orientation"]["Yz"];

	a[0][2] = a[1][0]*a[2][1] - a[2][0]*a[1][1];
	a[1][2] = a[2][0]*a[0][1] - a[0][0]*a[2][1];
	a[2][2] = a[0][0]*a[1][1] - a[1][0]*a[0][1];

	return apply_affine_transform(a,b,x);
}


int main (int argc, char ** argv)
{
    const int expected_argc = 2;
    if (argc < 1 + expected_argc)
    {
        std::cerr << "error: expected at least " <<  expected_argc << " command line arguments but detected " << argc << " of them:" << std::endl;
        for(auto i=0;i<argc; ++i)
            std::cerr << "argument #" << i << " = " << argv[i] << std::endl;
        return 1;
    }

	std::string metadata_filename(argv[1]);
	std::string query_point_filename(argv[2]);

	json metadata = read_metadata(metadata_filename);
	
	const uint64_t nx = metadata["width"];
    const uint64_t ny = metadata["height"];
	const uint64_t nz = metadata["depth"];
	assert(nx > 0);
	assert(ny > 0);
	assert(nz > 0);

	// create lattice of voxel-centered query points
	std::vector<std::pair<uint64_t, std::array<double, 3> > > query_points;
	for (uint64_t iz = 0; iz < nz; ++iz)
		for (uint64_t iy = 0; iy < ny; ++iy)
			for (uint64_t ix = 0; ix < nx; ++ix)
			{
				const uint64_t index = ix + nx*iy + nx*ny*iz;
				const std::array<double, 3> i({0.5+ix,0.5+iy,0.5+iz});
				const std::array<double, 3> transformed = apply_dicom_transform(metadata, i, false);
				std::pair<uint64_t, std::array<double, 3> > query_point;
				query_point.first = index;
				query_point.second = std::array<double, 3>({transformed[0], transformed[1], transformed[2]});
				query_points.push_back(query_point);
			}	

	// write to a text file
	std::ofstream outfile(query_point_filename);
	for (auto query_point : query_points) 
		outfile << query_point.first << " " << query_point.second[0] << " " << query_point.second[1] << " " << query_point.second[2] << std::endl;
	
	outfile.close();

    return 0;
}

