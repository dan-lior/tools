
#pragma once

#include "common_defs.h"
#include "affinity.h"

template<uint64_t dim>
struct Indexer
{
    Indexer(const Index<dim>& n);
    
    uint64_t to_rank(const Index<dim>& index) const;

    Index<dim> to_index(uint64_t rank) const;
    
    std::vector<uint64_t> neighbour_ranks(uint64_t rank, const Index<dim>& radii) const;
    
    bool out_of_range(Index<dim> index) const;
    
    uint64_t num_cells() const;
    
    const Index<dim> n;   
};

// Models a grid with arbitrary dimension, in a Euclidean space of arbitrary dimension
template<uint64_t dim_target, uint64_t dim_source>
struct Grid
{
    Grid(const Index<dim_source>& n, const Affinity<dim_target, dim_source>& affinity_) : indexer(n), affinity(affinity_) 
    {
    }
    
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

    const Indexer<dim_source> indexer;
    const Affinity<dim_target, dim_source> affinity;
};

// Models a grid with arbitrary dimension, embedded in a Euclidean space of the same dimension
template<uint64_t dim>
struct GridIso
{
    GridIso(const Index<dim>& n, const AffinityIso<dim>& affinity_ = AffinityIso<dim>()) : indexer(n), affinity(affinity_) 
    {
    }
    
    Vector<dim> position(const Vector<dim>& fractional_index) const
    {
        return affinity(fractional_index);
    }

    std::vector<Vector<dim>> all_positions() const
    {
        std::vector<Vector<dim>> rtn;
        for (uint64_t rank = 0; rank < indexer.num_cells(); rank++)
        {
            Index<dim> index = indexer.to_index(rank);
            Vector<dim> fractional_index = index_to_vector<dim>(index);
            Vector<dim> point = position(fractional_index);
            rtn.push_back(point);
        } 
        return rtn;
    }

    Vector<dim> fractional_index(const Vector<dim>& position) const
    {
        return affinity.inverse()(position);
    }

    // todo: compute fixed points!!!

    const Indexer<dim> indexer;
    const AffinityIso<dim> affinity;
};

// GridType is either Grid<dim_source, dim_target> or GridIso<dim>
template<typename GridType, typename T>
struct GridWithAttribute
{
    friend void export_to_vtk(const std::string& filename);
 
    GridWithAttribute(const GridType& grid_) :  grid(grid_)
    {
    }

    void set_attribute(const std::vector<T>& attribute_)
    {
        assert(attribute_.size() == grid.indexer.num_cells());
        attribute = attribute_;
    }

    std::vector<T> get_attribute() const
    {
        return attribute;
    }

    const GridType grid;

    private: 

    std::vector<T> attribute;   // don't ever change the size of this array!
};

namespace GridExport
{
    // the affinity stored in a grid is automatically applied to the points of the grid (via the position() member function),
    // the attribute, however, is NOT transformed. 

    template<typename GridType, typename T>
    void export_to_vtk(const GridWithAttribute<GridType, T>& gwa, const std::string& filename);
}


#include "voxel_box.hpp"