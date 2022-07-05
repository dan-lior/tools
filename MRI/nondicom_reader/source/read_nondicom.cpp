#include "read_nondicom.h"
#include "nlohmann/json.hpp"
#include "misc.h"
#include "affinity.h"
#include "voxel_box.h"

namespace
{
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

    template <typename T>
    double compute_speed_from_phase(const T& phase, double encoding_speed, uint64_t bit_depth)
    {
        // phase data ranges in [0,2^^bit_depth)
        // i ranges in [-2^^(bit_depth-1),2^^(bit_depth-1))
        // normalized phase data (n) ranges in [-1.0, 1.0)
        // velocity (v) components range in [-encoding_speed, encoding_speed)

        const uint64_t M = 1ULL<<(bit_depth-1);

        assert(phase <= 2*M);

        const int64_t i = static_cast<int64_t>(phase) - M;
        const double n = static_cast<double>(i) / M;
        const double v = encoding_speed*n;

        return v;
    }

    template <typename T>
    std::vector<Vector<3>> compute_velocities(const nlohmann::json& metadata, std::vector<T>& phasex, std::vector<T>& phasey, std::vector<T>& phasez)
    {
        AffinityIso<3> affinity = affinity_from_nondicom_metadata(metadata);
        Matrix<3> linear = affinity.linear_part;

        const uint64_t nx = metadata["width"];
        const uint64_t ny = metadata["height"];
        const uint64_t nz = metadata["depth"];

        const std::string temp(metadata["encoding_speed"]);
        const double encoding_speed = std::stod(temp);
        const uint64_t bit_depth = metadata["bit_depth"];

        std::vector<Vector<3>> rtn;

        for (uint64_t iz = 0; iz < nz; ++iz)
            for (uint64_t iy = 0; iy < ny; ++iy)
                for (uint64_t ix = 0; ix < nx; ++ix)
                {
                    uint64_t i = ny*nx*iz + nx*iy + ix;

                    Vector<3> v;
                    v(0) = compute_speed_from_phase(phasex[i], encoding_speed, bit_depth);
                    v(1) = compute_speed_from_phase(phasey[i], encoding_speed, bit_depth);
                    v(2) = compute_speed_from_phase(phasez[i], encoding_speed, bit_depth);

                    v = linear*v;
                    
                    v(2) *= -1; // hack! 

                    rtn.push_back(v);
                }	

        return rtn;
    }

    Index<3> grid_dimensions_from_nondicom_metadata(const nlohmann::json& nondicom_metadata)
    {
        const uint64_t nx = nondicom_metadata["width"];
        const uint64_t ny = nondicom_metadata["height"];
        const uint64_t nz = nondicom_metadata["depth"];
        assert(nx > 0);
        assert(ny > 0);
        assert(nz > 0);

        return Index<3>(nx,ny,nz);
    }

    GridIso<3> grid_from_nondicom_metadata(const nlohmann::json& nondicom_metadata)
    {
        return GridIso<3>(grid_dimensions_from_nondicom_metadata(nondicom_metadata), affinity_from_nondicom_metadata(nondicom_metadata));
    }

    std::vector<std::vector<Vector<3>> > get_velocity_from_phase_data_in_files(const nlohmann::json& metadata, const std::string& directory_name)
    {
        Index<3> grid_dimensions = grid_dimensions_from_nondicom_metadata(metadata);
        uint64_t num_cells = grid_dimensions(0) * grid_dimensions(1) * grid_dimensions(2);
        assert(num_cells > 0);

        // read phase data and write to vtk files
        assert(metadata["data_type"] == std::string("uint16")); // can only handle this datatype for now
        const std::string data_type = std::string("uint16");
        static_assert(sizeof(unsigned short)*8 == 16, "");
        typedef unsigned short T;

        const uint64_t num_timesteps = metadata["num_timesteps"];
        assert(num_timesteps > 0);

        std::vector<std::vector<T> > phasex_time_slices(num_timesteps);
        std::vector<std::vector<T> > phasey_time_slices(num_timesteps);
        std::vector<std::vector<T> > phasez_time_slices(num_timesteps);

        std::vector<std::vector<Vector<3>> > velocity_time_slices;

        for (uint64_t it=0; it<num_timesteps; ++it)
        {
            //debug
            std::cerr << "time step " << it << " of " << num_timesteps << std::endl;

            const std::string suffix(std::to_string(1+it));
            
            read_file(phasex_time_slices[it], directory_name + std::string("phasex_data/time_slices/") + suffix, num_cells);
            read_file(phasey_time_slices[it], directory_name + std::string("phasey_data/time_slices/") + suffix, num_cells);
            read_file(phasez_time_slices[it], directory_name + std::string("phasez_data/time_slices/") + suffix, num_cells);

            std::vector<Vector<3>> velocity_time_slice = compute_velocities(metadata, phasex_time_slices[it], phasey_time_slices[it], phasez_time_slices[it]);
            assert(velocity_time_slice.size() == num_cells);		

            velocity_time_slices.push_back(velocity_time_slice);
        }

        return velocity_time_slices;
    }
}

std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> read_nondicom(const std::string& nondicom_directory)
{
    const std::string metadata_filename = nondicom_directory + std::string("metadata.json");
    const nlohmann::json metadata = Miscellaneous::read_metadata(metadata_filename);

    GridIso<3> grid = grid_from_nondicom_metadata(metadata);

    std::vector<std::vector<Vector<3>> > velocity_timeslices = get_velocity_from_phase_data_in_files(metadata, nondicom_directory);
    assert(velocity_timeslices.size() == metadata["num_timesteps"]);

    std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> rtn;
    for (auto slice : velocity_timeslices)
    {
        GridWithAttribute<GridIso<3>, Vector<3>> temp(grid);
        temp.set_attribute(slice);
        rtn.push_back(temp);
    }
    return rtn;
}

