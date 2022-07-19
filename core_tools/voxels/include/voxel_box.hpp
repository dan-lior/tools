#include <cassert> 
#include "misc.h"

template<uint64_t dim_target, uint64_t dim_source>
Grid<dim_target, dim_source>::Grid(const Index<dim_source>& n, const Affinity<dim_target, dim_source>& affinity_) : indexer(n), affinity(affinity_) 
{}

template<uint64_t dim_target, uint64_t dim_source>
Vector<dim_target> Grid<dim_target, dim_source>::position(const Vector<dim_source>& fractional_index) const
{
    return affinity(fractional_index); 
}

template<uint64_t dim_target, uint64_t dim_source>
std::vector<Vector<dim_target>> Grid<dim_target, dim_source>::all_positions() const
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

template<uint64_t dim_target, uint64_t dim_source>
Vector<dim_source> Grid<dim_target, dim_source>::fractional_index(const Vector<dim_target>& position) const
{
    assert(affinity.isInvertible());
    return affinity.inverse()(position);
}

template<uint64_t dim_target, uint64_t dim_source, typename T>
LabelledGrid<dim_target, dim_source, T>::LabelledGrid(const Grid<dim_target, dim_source>& grid, const std::vector<T>& labels_) : 
    Grid<dim_target, dim_source>(grid), labels(labels_)
{
    using Base = Grid<dim_target, dim_source>;
    assert(labels.size() == Base::indexer.num_cells());
}

template<uint64_t dim_target, uint64_t dim_source, typename T>
std::vector<T> LabelledGrid<dim_target, dim_source, T>::get_labels() const
{
    return labels;
}

template<uint64_t dim_target, uint64_t dim_source, typename T>
LabelledGrid<1 + dim_target, 1 + dim_source, T> LabelledGrid<dim_target, dim_source, T>::stack_slices(const std::vector<LabelledGrid<dim_target, dim_source, T>>& slices, const Affinity<1,1>& time_map)
{
    const Grid<dim_target, dim_source> slice_grid = slices[0];

    MatrixRect<1 + dim_target, 1 + dim_source> A = MatrixRect<1 + dim_target, 1 + dim_source>::Zero();
    A.topLeftCorner(dim_target, dim_source) = slice_grid.affinity.linear_part;
    A.bottomRightCorner(1,1) = time_map.linear_part;

    Vector<1 + dim_target> b = Vector<1 + dim_target>::Zero();
    b(dim_target) = time_map.translation_part(0);
    
    Affinity<1 + dim_target, 1 + dim_source> affinity(A, b);

    Index<1 + dim_source> n = Index<1 + dim_source>::Zero();
    n.head(dim_source) = slice_grid.indexer.n;
    n(dim_source) = slices.size();

    Grid<1 + dim_target, 1 + dim_source> grid(n, affinity);

    std::vector<T> all_labels;
    for (auto slice : slices)
    {
        auto labels = slice.get_labels();
        all_labels.insert(all_labels.end(), labels.begin(), labels.end());
    }

    return LabelledGrid<1 + dim_target, 1 + dim_source, T>(grid, all_labels);
}


template<uint64_t dim_target, uint64_t dim_source, typename T>
std::vector<LabelledGrid<dim_target-1, dim_source-1, T>> LabelledGrid<dim_target, dim_source, T>::unstack_slices(const LabelledGrid<dim_target, dim_source, T>& stacked)
{
    // todo: handle the corner case: dim_target == 1 || dim_source == 1
    static_assert(dim_target > 1 && dim_source > 1, "slice too thin");

    MatrixRect<dim_target-1, dim_source-1> A = stacked.affinity.linear_part.topLeftCorner(dim_target-1, dim_source-1);
    Vector<dim_target-1> b = stacked.affinity.translation_part.head(dim_target-1);
    Affinity<dim_target-1, dim_source-1> affinity(A, b);
    Index<dim_source-1> n = stacked.indexer.n.head(dim_source-1);    
    Grid<dim_target-1, dim_source-1> slice_grid(n, affinity);

    const uint64_t num_slices = stacked.indexer.n(dim_source-1);
    assert(num_slices != 0 && stacked.indexer.num_cells() % num_slices == 0);
    const uint64_t slice_size = stacked.indexer.num_cells() / num_slices;

    const std::vector<T> labels = stacked.get_labels();

    std::vector<LabelledGrid<dim_target-1, dim_source-1, T>> rtn;
    for (uint64_t i=0; i<num_slices; ++i)
    {
        std::vector<T> slice_labels(labels.begin() + i*slice_size, labels.begin() + (i+1)*slice_size);
        rtn.push_back(LabelledGrid<dim_target-1, dim_source-1, T>(slice_grid, slice_labels));
    }

    return rtn;
}
