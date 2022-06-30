
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
    
    const Indexer<dim_source> indexer;
    const Affinity<dim_target, dim_source> affinity;
};

// Models a grid with arbitrary dimension, embedded in a Euclidean space of the same dimension
template<uint64_t dim>
struct GridIso
{
    GridIso(const Index<dim>& n, const AffinityIso<dim>& affinity_) : indexer(n), affinity(affinity_) 
    {
    }
    
    Vector<dim> position(const Vector<dim>& fractional_index) const
    {
        return affinity(fractional_index);
    }

    Vector<dim> fractional_index(const Vector<dim>& position) const
    {
        return affinity.inverse()(position);
    }

    // todo: compute fixed points!!!

    const Indexer<dim> indexer;
    const AffinityIso<dim> affinity;
};


Index<dim> grid_dimensions_from_dicom_metadata(const nlohmann::json& metadata)
{
    const uint64_t nx = metadata["width"];
    const uint64_t ny = metadata["height"];
    const uint64_t nz = metadata["depth"];

    return Index<dim>(nx,ny,nz);
}

AffinityIso<3> affinity_from_dicom_metadata(const nlohmann::json& metadata)
{
    // define a translation from the dicom data:
    Vector<3> b;
    b(0) = metadata["image_position"]["x"];
    b(1) = metadata["image_position"]["y"];
    b(2) = metadata["image_position"]["z"];

    // define a rotation from the dicom data:
    // according to the dicom spec, the first two columns are orthogonal and of unit length
    Matrix<3> a;

    a(0,0) = metadata["image_orientation"]["Xx"];
    a(1,0) = metadata["image_orientation"]["Xy"];
    a(2,0) = metadata["image_orientation"]["Xz"];

    a(0,1) = metadata["image_orientation"]["Yx"];
    a(1,1) = metadata["image_orientation"]["Yy"];
    a(2,1) = metadata["image_orientation"]["Yz"];

    a.col(2) = a.col(0).cross(a.col(1));

    // define a scaling from the dicom data
    const double dx = metadata["pixel_spacing"]["x"];
    const double dy = metadata["pixel_spacing"]["y"];
    const double dz = metadata["slice_thickness"];

    Matrix<3> d = Matrix<3>::Zero();
    d(0,0) = dx;
    d(1,1) = dy;
    d(2,2) = dz * (-1.0); // A hack to match conventions with 3d slicer

    return AffinityIso<3>(a*d, b);
}

GridIso<3> grid_from_dicom_metadata(const nlohmann::json& dicom_metadata)
{
    return GridIso(grid_dimensions_from_dicom_metadata(dicom_metadata), affinity_from_dicom_metadata(dicom_metadata));
}

// GridType is either Grid<dim_source, dim_target> or GridIso<dim>
template<typename GridType, typename T>
struct GridWithAttribute
{

    GridWithAttribute(const GridType& grid_) :  grid(grid_)
    {
    }

    void set_attribute(const std::vector<T>& attribute_)
    {
        assert(attribute_.size() == grid.indexer.num_cells());
        attribute = attribute_;
    }

    const GridType grid;
    std::vector<T> attribute;   // don't ever change the size of this array!
};


#include "voxel_box.hpp"