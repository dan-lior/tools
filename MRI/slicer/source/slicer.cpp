
#include "common_defs.h"
#include "helper.h"
#include "voxel_box.h"


using GWA = GridWithAttribute<GridIso<3>, Vector<3>>;

std::vector<GWA> read_voxel_box_time_slices(const nlohmann::json& metadata, const std::string& nondicom_directory)
{
	std::vector<std::vector<std::array<double, 3>> > velocity_timeslices = get_velocity_from_phase_data_in_files(metadata, nondicom_directory);
	assert(metadata["num_timesteps"] == velocity_timeslices.size());

	const AffinityIso<3> affinity = affinity_from_dicom_metadata(metadata); 
	const Index<3> n(metadata["height"], metadata["width"], metadata["depth"]);  //todo: check that this order is correct

	// transform velocities by the linear part of the dicom affinity
	std::vector<std::vector<Vector<3>> > transformed_velocity_timeslices;
	for (auto slice : velocity_timeslices)
	{
		std::vector<Vector<3>> transformed_slice;
		for (auto velocity : slice)
		{
			Vector<3> v(velocity[0], velocity[1], velocity[2]); // convert to a more convenient type
			Vector<3> tv = affinity.linear_part*v; // apply the linear part of the dicom affinity
			transformed_slice.push_back(tv);
		}

		transformed_velocity_timeslices.push_back(transformed_slice);
	}

	GWA volume_grid(GridIso<3>(n, affinity));
	std::vector<GWA> rtn;
	for (auto timeslice : transformed_velocity_timeslices)
	{
		volume_grid.set_attribute(timeslice);
		rtn.push_back(volume_grid);
	}

	return rtn;
}

template<typename T>
void sample_plane_in_volume(GridWithAttribute<Grid<3,2>, T>& plane_grid, const GridWithAttribute<GridIso<3>, T>& volume_grid)
{
	plane_grid.attribute.clear();

	for (uint64_t plane_rank=0; plane_rank < plane_grid.grid.indexer.num_cells(); ++plane_rank)
	{
		const Index<2> plane_index =  plane_grid.grid.indexer.to_index(plane_rank);
		const Vector<2> plane_fractional_index = index_to_vector<2>(plane_index);
		const Vector<3> pos = plane_grid.grid.position(plane_fractional_index);
		const Vector<3> volume_fractional_index = volume_grid.grid.fractional_index(pos);
		const Index<3> volume_index = vector_to_index<3>(volume_fractional_index);
		const uint64_t volume_rank = volume_grid.grid.indexer.to_rank(volume_index);
		const T atribute_value = volume_grid.attribute[volume_rank];

		plane_grid.attribute.push_back(atribute_value);
	}
}


int main (int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "error: expected 2 command line arguments but detected " << argc-1 << " of them:" << std::endl;
        for(auto i=1;i<argc; ++i)
            std::cerr << "command line argument #" << i << " = " << argv[i] << std::endl;
        return 1;
    }

	const std::string metadata_filename(argv[1]);
	const std::string nondicom_directory(argv[2]);
	nlohmann::json metadata = Miscellaneous::read_metadata(metadata_filename);


	// import velocity volume timeslices from MRI (appropriately transformed)
 	std::vector<GWA> volume_time_slices = read_voxel_box_time_slices(metadata, nondicom_directory);

	// specify sampling grid :

	// const double dx = 1.0;
	// const double dy = 1.0;
	// const double dz = 1.0;
	// const Vector<3> origin_of_plane;
	// const Vector<3> normal_to_plane;
	const Affinity<3,2> sampling_plane; // todo: specify this
	const uint64_t nx = 10;
	const uint64_t ny = 20;
	const Index<2> n({nx,ny});
	GridWithAttribute<Grid<3,2>, Vector<3>> sampling_grid(Grid<3,2>(n, sampling_plane));


	// sample volume along plane for each timeslice

	std::vector<GridWithAttribute<Grid<3,2>, Vector<3>>> planar_timeslices;
	for (auto volume_timeslice : volume_time_slices)
	{		
		sample_plane_in_volume(sampling_grid, volume_timeslice);
		planar_timeslices.push_back(sampling_grid);
	}

    return 0;
}
