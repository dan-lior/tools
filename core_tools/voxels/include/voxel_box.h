
#pragma once

#include "common_defs.h"
#include "affinity.h"
#include "indexer.h"

template<uint64_t dim_target, uint64_t dim_source>
struct Grid
{
    Grid(const Index<dim_source>& n, const Affinity<dim_target, dim_source>& affinity_) : indexer(n), affinity(affinity_) 
    {}
    
    Vector<dim_target> position(const Vector<dim_source>& fractional_index) const
    {
       return affinity(fractional_index); 
    }
    
    std::vector<Vector<dim_target>> all_positions() const
    {
        std::vector<Vector<dim_target>> rtn;
        for (uint64_t rank = 0; rank < indexer.num_cells(); rank++)
        {
            Index<dim_source> index = indexer.to_index(rank);
            Vector<dim_source> fractional_index = index_to_vector<dim_source>(index);
            Vector<dim_target> point = position(fractional_index);
            rtn.push_back(point);
        } 
        return rtn;
    }

    Vector<dim_source> fractional_index(const Vector<dim_target>& position) const
    {
        assert(affinity.isInvertible());

        return affinity.inverse()(position);
    }

    const Indexer<dim_source> indexer;
    const Affinity<dim_target, dim_source> affinity;
};

template<uint64_t dim_target, uint64_t dim_source, typename T>
struct LabelledGrid : public Grid<dim_target, dim_source>
{
    LabelledGrid(const Grid<dim_target, dim_source>& grid, const std::vector<T>& labels_) : 
        Grid<dim_target, dim_source>(grid), labels(labels_)
    {}

    // during export, the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the attribute itself, however, is NOT transformed. 
    void export_to_vtk(const std::string& filename);

    std::vector<T> get_labels() const
    {
        return labels;
    }

    private: 
    const std::vector<T> labels;
};
