
#include "disjoint_union.h"

DisjointUnion::DisjointUnion(uint64_t n)
{
    parents.resize(n);
    for (uint64_t i=0; i<n; ++i)
        parents[i]=i;
}

uint64_t DisjointUnion::number_of_disjoint_subsets() const
{
    uint64_t num_components = 0;
    for (uint64_t i=0; i<parents.size(); ++i)
        if (i == get_representative(i))
            num_components++;

    return num_components;
}

std::vector<uint64_t> DisjointUnion::disjoint_subset(uint64_t i) const
{
    assert(i < parents.size());
    uint64_t ii = get_representative(i);

    std::vector<uint64_t> rtn;
    for (uint64_t j=0; j<parents.size(); ++j)
        if (ii == get_representative(j))
            rtn.push_back(j);

    return rtn;            
}

uint64_t DisjointUnion::get_representative(uint64_t i) const 
{
    assert(i<parents.size());
    while(i != parents[i])
        i = parents[i];

    return i;
}

void DisjointUnion::merge(uint64_t i, uint64_t j)
{
    assert(i<parents.size());
    assert(j<parents.size());
    uint64_t ii = get_representative(i);
    uint64_t jj = get_representative(j);
    if (ii < jj)
        parents[ii] = jj;
    else if (ii > jj)
        parents[jj] = ii;
}

void DisjointUnion::cleanup() // invoking this, with the appropriate frequency, makes the data structure more efficient
{
    for (uint64_t i=0; i<parents.size(); ++i)
        parents[i] = get_representative(i);
}

