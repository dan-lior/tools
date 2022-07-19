#pragma once

#include "common_defs.h"

template<uint64_t dim>
struct Indexer
{
    Indexer(const Index<dim>& n);
    
    uint64_t to_rank(const Index<dim>& index) const;

    Index<dim> to_index(uint64_t rank) const;
    
    std::vector<uint64_t> neighbour_ranks(uint64_t rank, const Index<dim>& radii) const;
    
    bool out_of_range(Index<dim> index) const;
    
    uint64_t num_cells() const;
    
    Index<dim> n;   
};

#include "indexer.hpp"