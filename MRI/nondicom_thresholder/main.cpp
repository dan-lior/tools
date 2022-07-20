
#include <iostream>
#include <vector>
#include <cassert>
#include "misc.h"
#include "voxel_box.h"
#include "read_nondicom.h"

namespace
{
    std::vector<double> all_speeds(const std::vector<LabelledGrid<3, 3, Vector<3>> >& velocity_time_slices)
    {
        std::vector<double> rtn;

        for (auto velocity_time_slice : velocity_time_slices)
        {
            auto velocities = velocity_time_slice.get_labels();
            for (auto velocity : velocities)
                rtn.push_back(velocity.norm());

            std::cerr << "num speeds: " << rtn.size();
            std::cerr << " min speed: " << *std::min_element(rtn.begin(), rtn.end());
            std::cerr << " max speed: " << *std::max_element(rtn.begin(), rtn.end()) << std::endl;
        }

        return rtn;
    }

    std::vector<double> all_mean_cosines(const std::vector<LabelledGrid<3, 3, Vector<3>> >& velocity_time_slices, 
    uint64_t rx, uint64_t ry, uint64_t rz, uint64_t rt)
    {
        const uint64_t nt = velocity_time_slices.size();

        const Index<3> n = velocity_time_slices[0].indexer.n;
        const uint64_t nx = n(0);
        const uint64_t ny = n(1);
        const uint64_t nz = n(2);

        assert(nt > rt);
        assert(nx > rx);
        assert(ny > ry);
        assert(nz > rz);

        std::vector<double> speeds = all_speeds(velocity_time_slices);

        std::vector<double> mean_cosines;

        for (uint64_t it = 0; it < nt; ++it)
        {
            auto labels = velocity_time_slices[it].get_labels();

            std::cerr << "computing mean cosines for slice: " << it << std::endl;
            for (uint64_t iz = 0; iz < nz; ++iz)
                for (uint64_t iy = 0; iy < ny; ++iy)
                    for (uint64_t ix = 0; ix < nx; ++ix)
                    {

                        std::vector<double> cos_thetas;
                        const uint64_t index = ix + nx*iy + nx*ny*iz;
                        const uint64_t big_index = index + it*nx*ny*nz;

                        // loop over neighbours
                        uint64_t min_jt = it < rt  ? 0 : it-rt;
                        uint64_t max_jt = std::min<uint64_t>(nt, it+rt+1);
                        assert(max_jt - min_jt >= rt);


                        for (uint64_t jt = min_jt; jt < max_jt; ++jt)
                        {
                            uint64_t min_jz = iz < rz ? 0 : iz -rz;
                            uint64_t max_jz = std::min<uint64_t>(nz, iz+rz+1);
                            assert(max_jz - min_jz >= rz);

                            for (uint64_t jz = min_jz; jz < max_jz; ++jz)
                            {
                                uint64_t min_jy = iy < ry ? 0 : iy-ry;
                                uint64_t max_jy = std::min<uint64_t>(ny, iy+ry+1);
                                assert(max_jy - min_jy >= ry);

                                for (uint64_t jy = min_jy; jy < max_jy; ++jy)
                                {
                                    uint64_t min_jx = ix < rx ? 0 : ix -rx;
                                    uint64_t max_jx = std::min<uint64_t>(nx, ix+rx+1);
                                    assert(max_jx - min_jx >= rx);

                                    for (uint64_t jx = min_jx; jx < max_jx; ++jx)
                                    {
                                        const uint64_t index_neighbour = jx + nx*jy + nx*ny*jz;
                                        const uint64_t big_index_neighbour = index_neighbour + jt*nx*ny*nz;

                                        const double dot = (labels[index].dot(labels[index_neighbour]));

                                        double cos_theta = 0;
                                        if (dot != 0)
                                        {
                                            auto denom = speeds[big_index] * speeds[big_index_neighbour];
                                            assert(denom != 0);
                                            cos_theta = dot / denom;
                                        }

                                        cos_thetas.push_back(cos_theta);
                                    }
                                }
                            }
                        }

                        assert(cos_thetas.size() >= (1+rx) * (1+ry) * (1+ rz)* (1+ rt));
                        assert(cos_thetas.size() <= (1+2*rx) * (1+2*ry) * (1+ 2*rz)* (1+ 2*rt));

                        double sum_cosine = 0;
                        for (auto cosine : cos_thetas)
                            sum_cosine += cosine;
                        mean_cosines.push_back(sum_cosine / cos_thetas.size());

                    }	
        }

        assert(mean_cosines.size() == nx*ny*nz*nt);
        return mean_cosines;

    }

    /*
    std::vector<double> all_speed_variances(const nlohmann::json& metadata, const std::vector<std::vector<std::array<double, 3>> >& velocity_time_slices, uint64_t rx, uint64_t ry, uint64_t rz, uint64_t rt)
    {
        std::vector<uint64_t> rtn;

        const uint64_t nx = metadata["width"];
        const uint64_t ny = metadata["height"];
        const uint64_t nz = metadata["depth"];
        const uint64_t nt = metadata["num_timesteps"];

        assert(nt > rt);
        assert(nx > rx);
        assert(ny > ry);
        assert(nz > rz);

        assert(velocity_time_slices.size() == nt);
        for (auto slice : velocity_time_slices)
            assert(slice.size() == nx*ny*nz);

        std::vector<double> speeds = all_speeds(velocity_time_slices);

        // for each voxel (in time and space) compute the variance of the speed in a neighbourhood around that voxel
        std::vector<double> all_variances;

        for (uint64_t it = 0; it < nt; ++it)
        {
            std::cerr << "computing speed variances for slice: " << it << std::endl;
            for (uint64_t iz = 0; iz < nz; ++iz)
                for (uint64_t iy = 0; iy < ny; ++iy)
                    for (uint64_t ix = 0; ix < nx; ++ix)
                    {
                        std::vector<double> neighbourhood_speeds;

                        // loop over neighbours
                        uint64_t min_jt = it < rt  ? 0 : it-rt;
                        uint64_t max_jt = std::min<uint64_t>(nt, it+rt+1);
                        assert(max_jt - min_jt >= rt);
                        
                        for (uint64_t jt = min_jt; jt < max_jt; ++jt)
                        {
                            uint64_t min_jz = iz < rz ? 0 : iz -rz;
                            uint64_t max_jz = std::min<uint64_t>(nz, iz+rz+1);
                            assert(max_jz - min_jz >= rz);

                            for (uint64_t jz = min_jz; jz < max_jz; ++jz)
                            {
                                uint64_t min_jy = iy < ry ? 0 : iy-ry;
                                uint64_t max_jy = std::min<uint64_t>(ny, iy+ry+1);
                                assert(max_jy - min_jy >= ry);

                                for (uint64_t jy = min_jy; jy < max_jy; ++jy)
                                {
                                    uint64_t min_jx = ix < rx ? 0 : ix -rx;
                                    uint64_t max_jx = std::min<uint64_t>(nx, ix+rx+1);
                                    assert(max_jx - min_jx >= rx);

                                    for (uint64_t jx = min_jx; jx < max_jx; ++jx)
                                    {
                                        const uint64_t index_neighbour = jx + nx*jy + nx*ny*jz + nx*ny*nz*jt;
                                        neighbourhood_speeds.push_back(speeds[index_neighbour]);
                                    }
                                }
                            }
                        }

                        all_variances.push_back(variance(neighbourhood_speeds));
                    }	
        }

        assert(all_variances.size() == nx*ny*nz*nt);
        return all_variances;
    }

    */

    std::vector<uint64_t> threshold(const std::vector<double>& values, double threshold, bool above_not_below)
    {
        std::cerr << "min value: " << *std::min_element(values.begin(), values.end()) << std::endl;
        std::cerr << "max value: " << *std::max_element(values.begin(), values.end()) << std::endl;
        std::cerr << "threshold: " << threshold << std::endl;

        std::vector<uint64_t> selected_indices;

        for (uint64_t i=0; i<values.size(); ++i)
            if ((values[i] > threshold) == above_not_below)
                selected_indices.push_back(i);

        std::cerr << "percent within threshold: " << 1.0*selected_indices.size() / values.size() << std::endl;

        return selected_indices;
    }


    std::vector<uint64_t> compute_outside_indices(const std::string& mask_index_filename, uint64_t num_voxels_per_timeslice,  uint64_t num_timeslices)
    {
        // read inside indices, take complement, and repeat over time slices
        std::vector<uint64_t> outside_indices;
        std::set<uint64_t> inside_indices;

        std::ifstream infile(mask_index_filename);
        std::string str;
        while(std::getline(infile, str))
        {
            uint64_t index = std::stoi(str);
            inside_indices.insert(index);
        }
        infile.close();
        std::cerr << "number of inside indices in slice: " << inside_indices.size() << std::endl;

        std::vector<uint64_t> outside_indices_slice;
        for (uint64_t index = 0; index < num_voxels_per_timeslice; ++index)
            if (inside_indices.find(index) == inside_indices.end())
                outside_indices_slice.push_back(index);
        std::cerr << "num outside indices in slice: " << outside_indices_slice.size() << std::endl;
            
        assert(inside_indices.size() + outside_indices_slice.size() == num_voxels_per_timeslice);

        for (uint64_t it=0; it<num_timeslices; ++it)
            for (uint64_t i : outside_indices_slice)
                outside_indices.push_back(it*num_voxels_per_timeslice+ i);
        std::cerr << "num outside indices : " << outside_indices.size() << std::endl;

        assert(outside_indices.size() == num_timeslices * outside_indices_slice.size());

        return outside_indices;
    }

    template<typename T>
    std::vector<T> clear_labels(std::vector<T> labels, const std::vector<uint64_t>& ranks)
    {
        for (auto rank : ranks)
            labels[rank] = T();

        return labels;
    }
}


int main()
{
    const std::string nondicom_directory("/home/dan/git-repos/tools/MRI/nondicom_data/patient1/");
    assert(Miscellaneous::directory_exists(nondicom_directory));

	const std::string mask_index_filename("/home/dan/git-repos/tools/MRI/nondicom_thresholder/data/blah.txt");
    assert(Miscellaneous::file_exists(mask_index_filename));

	const std::string vtk_directory("/home/dan/git-repos/tools/MRI/nondicom_thresholder/vtk/");
    std::filesystem::create_directory(vtk_directory);

	const std::string vtk_test1_directory = vtk_directory + std::string("test1/");
    std::filesystem::create_directory(vtk_test1_directory);

	const std::string vtk_test2_directory = vtk_directory + std::string("test2/");
    std::filesystem::create_directory(vtk_test2_directory);

	const std::string vtk_mask_directory = vtk_directory + std::string("masks/");
    std::filesystem::create_directory(vtk_mask_directory);

	const std::string vtk_velocity_directory = vtk_directory + std::string("velocities/");
    std::filesystem::create_directory(vtk_velocity_directory);

    std::vector<LabelledGrid<3, 3, Vector<3>>> velocity_timeslices_all = read_nondicom(nondicom_directory);
    

    std::vector<LabelledGrid<3, 3, Vector<3>>> velocity_timeslices;
    velocity_timeslices.push_back(velocity_timeslices_all[0]);
    velocity_timeslices.push_back(velocity_timeslices_all[1]);
    velocity_timeslices.push_back(velocity_timeslices_all[2]);
    velocity_timeslices.push_back(velocity_timeslices_all[3]);

    // for (uint64_t i=0; i<velocity_timeslices.size(); ++i)
    //     LabelledGrid<3, 3, Vector<3>> :: export_labelled_grid_to_vtk(velocity_timeslices[i], vtk_test1_directory + std::to_string(i) + std::string(".vtk"));

    MatrixRect<1,1> A;
    A(0,0) = 37.5; // time_between_slices, in milliseconds;
    Vector<1> b;
    b(0) = 88.8; // time_at_first_slice;

    Affinity<1,1> time_map(A,b); 
    LabelledGrid<4, 4, Vector<3>> stacked = LabelledGrid<3, 3, Vector<3>>::stack_slices(velocity_timeslices, time_map);

    // std::vector<LabelledGrid<3, 3, Vector<3>>> unstacked = LabelledGrid<4, 4, Vector<3>>::unstack_slices(stacked);

    // for (uint64_t i=0; i<unstacked.size(); ++i)
    //     LabelledGrid<3, 3, Vector<3>> :: export_labelled_grid_to_vtk(unstacked[i], vtk_test2_directory + std::to_string(i) + std::string(".vtk"));

    const Index<4> nhd_dimensions(2,2,2,1);   //rx, ry, rz, rt
    std::vector<std::vector<Index<4>>> neighbourhoods = stacked.indexer.neighbourhoods(nhd_dimensions);
    
    std::cerr << "nhd info:" << std::endl;
    for (auto nhd : neighbourhoods)
    {
        std::cerr << nhd.size() << " ";
    }
    std::cerr << std::endl;



    return 0; //debug


    std::vector<double> mean_cosines = all_mean_cosines(velocity_timeslices, 2,2,2,1); //rx, ry, rz, rt



	const uint64_t nt = velocity_timeslices.size();
    const uint64_t n = velocity_timeslices[0].indexer.num_cells();


    const double low_mean_cosine_threshold = 0.7;
	std::vector<uint64_t> indices_beyond_threshold_mean_cosine = threshold(mean_cosines, low_mean_cosine_threshold, false);

	std::vector<uint64_t> outside_indices = compute_outside_indices(mask_index_filename, n, nt); // indices of voxels outside the segmented surface

	std::set<uint64_t> united_indices;
	united_indices.insert(indices_beyond_threshold_mean_cosine.begin(), indices_beyond_threshold_mean_cosine.end());
	united_indices.insert(outside_indices.begin(), outside_indices.end());

    std::vector<std::vector<uint64_t>> indices_slice_by_slice(nt);
    for (auto index : united_indices)
    {
        uint64_t it = index / n;
        uint64_t index_within_slice = index % n;
        indices_slice_by_slice[it].push_back(index_within_slice);
    }

    std::vector<LabelledGrid<3, 3, uint64_t >> mask_timeslices;

    const Grid<3, 3> grid = velocity_timeslices[0];

    for (uint64_t it=0; it<velocity_timeslices.size(); ++it)
    {
        // const Grid<3, 3> grid = velocity_timeslices[it]; // in fact, this grid is the same for all it

        auto cleared_velocity_labels = clear_labels(velocity_timeslices[it].get_labels(), indices_slice_by_slice[it]);
        velocity_timeslices[it] = LabelledGrid<3, 3, Vector<3>>(grid, cleared_velocity_labels);
        LabelledGrid<3, 3, Vector<3>>::export_labelled_grid_to_vtk(velocity_timeslices[it], vtk_velocity_directory + std::to_string(it) + std::string(".vtk"));

        auto cleared_mask_labels = clear_labels(std::vector<uint64_t> (n, 1), indices_slice_by_slice[it]);
        mask_timeslices.push_back(LabelledGrid<3, 3, uint64_t>(grid, cleared_mask_labels));
        LabelledGrid<3, 3, uint64_t>::export_labelled_grid_to_vtk(mask_timeslices[it], vtk_mask_directory + std::to_string(it) + std::string(".vtk"));
    }

    return 0;
}

