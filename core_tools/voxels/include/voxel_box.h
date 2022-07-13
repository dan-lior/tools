
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

    const Indexer<dim_source> indexer;
    const Affinity<dim_target, dim_source> affinity;
};

template<uint64_t dim_target, uint64_t dim_source, typename T>
struct LabelledGrid : public Grid<dim_target, dim_source>
{
    LabelledGrid(const Grid<dim_target, dim_source>& grid, const std::vector<T>& labels);

    std::vector<T> get_labels() const;

    // during export, the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the attribute itself, however, is NOT transformed. 
    void export_to_vtk(const std::string& filename) const;

    private: 

    const std::vector<T> labels;

};


#include "voxel_box.hpp"