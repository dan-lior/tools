
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

    const Indexer<dim_source> indexer;
    const Affinity<dim_target, dim_source> affinity;
};

template<uint64_t dim_target, uint64_t dim_source, typename T>
struct LabelledGrid : public Grid<dim_target, dim_source>
{
    LabelledGrid(const Grid<dim_target, dim_source>& grid, const std::vector<T>& labels);

    std::vector<T> get_labels() const;

    // query this object to generate labels for probe
    template<uint64_t dim_source_probe>
    LabelledGrid<dim_target, dim_source_probe, T> label_probe(const Grid<dim_target, dim_source_probe>& probe) const;

    // return a copy of the curent object with all labels cleared except those specified by probe
    template<uint64_t dim_source_probe>
    LabelledGrid<dim_target, dim_source, T> clear_labels_outside_probe(const Grid<dim_target, dim_source_probe>& probe) const;

    private: 

    const std::vector<T> labels;
};


#include "voxel_box.hpp"