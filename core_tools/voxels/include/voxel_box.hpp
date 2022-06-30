#include "misc_math.h"


template<uint64_t dim>
Indexer<dim>::Indexer(const Index<dim>& n_) : n(n_) {}

template<uint64_t dim>
uint64_t Indexer<dim>::to_rank(const Index<dim>& index) const
{
    uint64_t prod = 1;
    uint64_t rank = 0;
    for (uint64_t i=0; i<dim; ++i)
    {
        rank += prod*index(i);
        prod *= n(i);
    }
    return rank;
}

template<uint64_t dim>
Index<dim> Indexer<dim>::to_index(uint64_t rank) const
{
    Index<dim> index = Index<dim>::Zero();
    for (uint64_t i=0; i<dim; ++i)
    {
        index(i) = rank % n(i);
        rank = (rank -index(i))/n(i);
    }
    return index;
}

template<uint64_t dim>
std::vector<uint64_t> Indexer<dim>::neighbour_ranks(uint64_t rank, const Index<dim>& radii) const
{
    assert(rank < num_cells);
    Index<dim> index = to_index(rank);

    std::vector<std::vector<uint64_t>> intervals;
    for (uint64_t d=0; d<dim; ++d)
    {
        const uint64_t r = radii[d];
        std::vector<uint64_t> interval;
        int64_t i = index(d);

        // uint64_t i_beg = std::max<uint64_t>(0,i-r);
        uint64_t i_beg = i>r ? i-r : 0;
        uint64_t i_end = std::min<uint64_t>(i+r+1, n(d));
        for (int64_t i=i_beg; i<i_end; ++i)
            interval.push_back(i);
        intervals.push_back(interval);
    }

    std::vector<uint64_t> rtn;

    auto indices = MiscellaneousMath::cartesian_product(intervals);
    for (auto index : indices)
    {
        assert(index.size() == dim);
        Index<dim> index_;
        for (uint64_t i=0; i<dim; ++i)
            index_(i) = index[i];
        uint64_t rank = to_rank(index_);
        assert(rank < num_cells);
        rtn.push_back(rank);
    }
    
    return rtn;
}

template<uint64_t dim>
bool Indexer<dim>::out_of_range(Index<dim> index) const
{
    for (uint64_t i=0; i<dim; ++i)
        if (index(i) < 0 || index(i) >= n(i))
            return true;

    return false;
}

template<uint64_t dim>
uint64_t Indexer<dim>::num_cells() const
{
    uint64_t rtn = 1;
    for (uint64_t i=0; i<dim; ++i)
        rtn *= n(i);

    return rtn;
}
