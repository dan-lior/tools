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
void LabelledGrid<dim_target, dim_source, T>::export_to_vtk(const std::string& filename) const
{
    std::cerr << "export_labelled_grid_to_vtk is only implemented for certain template specializations ... not the ones supplied " << std::endl;
}

