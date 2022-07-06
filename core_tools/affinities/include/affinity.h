#pragma once

#include "common_defs.h"


// models an affine map: x --> Ax + b

template<uint64_t dim_target, uint64_t dim_source>
struct Affinity
{
    static constexpr uint64_t target_dimension = dim_target;
    static constexpr uint64_t source_dimension = dim_source;

    Affinity(
        const MatrixRect<dim_target, dim_source>& linear = MatrixRect<dim_target, dim_source>::Identity(), 
        const Vector<dim_target>& translation = Vector<dim_target>::Zero());
    
    // eg, for dim_source=2, this matrix acts on a vector of the form [x y 1]^t
    MatrixRect<1+dim_target, 1+dim_source> single_matrix_representation() const;

    // composition (A*B)(x) := A(B(x))
    template<uint64_t dim_presource>
    Affinity<dim_target, dim_presource> operator*(const Affinity<dim_target, dim_source>& other) const;

    // evaluation
    Vector<dim_target> operator()(const Vector<dim_source>& p) const;

    bool isInvertible() const;

    Affinity<dim_source, dim_target> inverse() const;

    static Affinity<dim_target, dim_source> procrustes(const std::vector<Vector<dim_source>>& source, const std::vector<Vector<dim_target>>& target, bool force_special_orthogonal = false)
    {
        if (dim_target == dim_source)
            std::cerr << "not yet implemented" << std::endl;
        else
            std::cerr << "not defined" << std::endl;

        return Affinity<dim_target, dim_source>();
    }

    // return false if and only if no fixed points exist
    // the set of fixed points have the form particular_fixed_point + linear combo of basis elements
    bool compute_fixed_points(Vector<dim_target>& particular_fixed_point, std::vector<Vector<dim_target>>& basis) const
    {
        assert(dim_target == dim_source);

        std::cerr << "not yet implemented" << std::endl;

        return true;
    }

    MatrixRect<dim_target, dim_source> linear_part;
    Vector<dim_target> translation_part;
};


// temporary hack ... remove soon
template<uint64_t dim>
using AffinityIso = Affinity<dim, dim>;

#include "affinity.hpp"