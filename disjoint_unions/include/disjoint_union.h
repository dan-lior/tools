
#pragma once

#include "common_defs.h"

struct DisjointUnion
{
    DisjointUnion(uint64_t n);

    uint64_t number_of_disjoint_subsets() const;

    std::vector<uint64_t> disjoint_subset(uint64_t i) const;

    uint64_t get_representative(uint64_t i) const;
    
    void merge(uint64_t i, uint64_t j);

    void cleanup(); // invoking this, with the appropriate frequency, makes the data structure more efficient

    std::vector<uint64_t> parents;
};
