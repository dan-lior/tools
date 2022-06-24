#pragma once

#include "voxel_box.h"



template<uint64_t dim>
class VoxelModel
{
    public: 

    VoxelModel<dim>(const VoxelBox<dim>& voxel_box_, const std::vector<bool>& voxel_mask_);
    
    static void display_mask_stats(const std::vector<bool>& mask);

    // VoxelModel(const VoxelBox& voxel_box, const std::vector<uint64_t>& voxel_ranks);
    // VoxelModel(const VoxelBox& voxel_box, const std::vector<std::array<uint64_t, 3> >& voxel_indices);

    // Index<dim> neighbourhood; // specify the radii of a box
    // VoxelModel get_interior(const Neighbourhood& nhd) const;
    // VoxelModel get_boundary(const Neighbourhood& nhd) const;

    std::vector<bool> transformed_mask(const Affinity<dim> affinity) const;


    std::vector<bool> boundary_mask(const Index<dim>& radii) const;
    std::vector<bool> plane_cut_mask(const Vector<dim>& point, const Vector<dim>& normal);

    // returns the index of a closest voxel in other to the voxel with the specified index
    Index<dim> closest_point(Index<dim> index, const VoxelModel<dim>& other) const;

    void write_to_vtk(const std::string& filename) const;

    std::vector<bool> voxel_mask; 
    const VoxelBox<dim> voxel_box;
};

#include "voxel_model.hpp"