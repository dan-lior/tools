
#pragma once

#include "common_defs.h"
#include "affinity.h"
#include "indexer.h"



template<uint64_t dim_target, uint64_t dim_source>
struct Grid
{
    Grid(const Index<dim_source>& n, const Affinity<dim_target, dim_source>& affinity);
    
    Vector<dim_target> position(const Vector<dim_source>& fractional_index) const;
    
    std::vector<Vector<dim_target>> all_positions() const;

    Vector<dim_source> fractional_index(const Vector<dim_target>& position) const;

    std::vector<std::vector<uint64_t>> neighbourhoods(const Index<dim_source>& nhd_dimensions) const
    {
        std::vector<std::vector<uint64_t>> rtn;
        for (uint64_t rank=0; rank<indexer.num_cells(); ++rank)
        {
            Index<dim_source> index = indexer.to_index(rank);

            Index<dim_source> min_corner;
            Index<dim_source> max_corner;
            for (uint64_t i = 0; i < dim_source; ++ i)
            {
                max_corner(i) = std::max(index(i) + nhd_dimensions[i], (indexer.n)(i));
                min_corner(i) = (i >= nhd_dimensions[i]) ? i - nhd_dimensions[i] : 0;
                std::max(i + nhd_dimensions[i], (indexer.n)(i));
            }

            std::vector<uint64_t> nhd_of_indices = MiscellaneousMath::cartesian_product(min_corner, max_corner);



            std::vector<uint64_t> nhd_of_ranks;

            rtn.push_back(nhd_of_ranks);
        }
        return rtn;
    }





    // other grid must have an invertible affinity
    // returns an associative map (generally not 1-1 nor onto)
    // (r,s) is in the map if center of voxel r is inside voxel s
    
    std::map<uint64_t, uint64_t> rank_map(const Grid<dim_target, dim_target>& other) const
    {
        assert(other.affinity.isInvertible());

        std::map<uint64_t, uint64_t> rtn;
        for (uint64_t rank = 0; rank < indexer.num_cells(); ++rank)
        {
            Index<dim_source> index = indexer.to_index(rank);
            Vector<dim_source> fractional_index = index_to_vector<dim_source>(index) + 0.5 * Vector<dim_source>::Ones();
            Vector<dim_target> pos = position(fractional_index);
            Vector<dim_target> fractional_index_other = other.fractional_index(pos); // this step needs invertibility
            Index<dim_target> index_other = vector_to_index<dim_target>(fractional_index_other);

            if (!other.indexer.out_of_range(index_other))
                rtn[rank] = other.indexer.to_rank(index_other);
        }
        return rtn;
    }

    Indexer<dim_source> indexer;
    Affinity<dim_target, dim_source> affinity;
};

template<uint64_t dim_target, uint64_t dim_source, typename T>
struct LabelledGrid : public Grid<dim_target, dim_source>
{
    LabelledGrid(const Grid<dim_target, dim_source>& grid, const std::vector<T>& labels);
    
    std::vector<T> get_labels() const;


    // NOTE: slices are always in the last dimension (ie for an xyz cube, slices are parallel to the xy plane)

    // assumes that all the slices have the same underlying grid
    // time map specifies location and spacing of slices
    static LabelledGrid<1 + dim_target, 1 + dim_source, T> stack_slices(const std::vector<LabelledGrid<dim_target, dim_source, T>>& slices, const Affinity<1,1>& time_map = Affinity<1,1>());

    // discards information about time map
    static std::vector<LabelledGrid<dim_target-1, dim_source-1, T>> unstack_slices(const LabelledGrid<dim_target, dim_source, T>& stacked_slices);


    // during export, the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the labels themselves, however, are NOT transformed. 
    static void export_labelled_grid_to_vtk(const LabelledGrid<dim_target, dim_source, T>& labelled_grid, const std::string& filename);

private: 

    std::vector<T> labels;
};


#include "voxel_box.hpp"