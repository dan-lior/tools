#include <cassert> 

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
template<uint64_t dim_source_probe>
LabelledGrid<dim_target, dim_source_probe, T> LabelledGrid<dim_target, dim_source, T>::label_probe(const Grid<dim_target, dim_source_probe>& probe) const
{
    using Base = Grid<dim_target, dim_source>;

    assert(Base::affinity.isInvertible()); 
    // in particular dim_source == dim_target

    std::vector<T> probe_labels(probe.indexer.num_cells());

    for (uint64_t probe_rank = 0; probe_rank < probe.indexer.num_cells(); ++probe_rank)
    {
        const Index<dim_source_probe> probe_index = probe.indexer.to_index(probe_rank);
        const Vector<dim_source_probe> probe_fractional_index = index_to_vector<dim_source_probe>(probe_index);
        const Vector<dim_target> pos = probe.position(probe_fractional_index);
        const Vector<dim_source> fractional_index = Base::fractional_index(pos); // this is where invertibility is needed
        const Index<dim_source> index = vector_to_index<dim_source>(fractional_index);

        if (!Base::indexer.out_of_range(index))
            probe_labels[probe_rank] = labels[Base::indexer.to_rank(index)];
    }

    return LabelledGrid<dim_target, dim_source_probe, T>(probe, probe_labels);
}


template<uint64_t dim_target, uint64_t dim_source, typename T>
template<uint64_t dim_source_probe>
LabelledGrid<dim_target, dim_source, T> LabelledGrid<dim_target, dim_source, T>::clear_labels_outside_probe(const Grid<dim_target, dim_source_probe>& probe) const
{
    using Base = Grid<dim_target, dim_source>;

    assert(Base::affinity.isInvertible()); 
    // in particular dim_source == dim_target

    std::vector<T> new_labels(labels.size(), Vector<3>());
    
    for (uint64_t probe_rank = 0; probe_rank < probe.indexer.num_cells(); ++probe_rank)
    {
        const Index<dim_source_probe> probe_index = probe.indexer.to_index(probe_rank);
        const Vector<dim_source_probe> probe_fractional_index = index_to_vector<dim_source_probe>(probe_index);
        const Vector<dim_target> pos = probe.position(probe_fractional_index);
        const Vector<dim_source> fractional_index = Base::fractional_index(pos); // this is where invertibility is needed
        const Index<dim_source> index = vector_to_index<dim_source>(fractional_index);

        if (!Base::indexer.out_of_range(index))
            new_labels[probe_rank] = labels[Base::indexer.to_rank(index)];
    }

    return LabelledGrid<dim_target, dim_source, T>(Base(Base::indexer.n, Base::affinity), new_labels);
}




