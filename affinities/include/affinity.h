#pragma once

#include "common_defs.h"



template<uint64_t dim_target, uint64_t dim_source>
struct Affinity
{
    // specifies the affine map: x --> Ax + b
    Affinity(
        const MatrixRect<dim_target, dim_source>& linear = MatrixRect<dim_target, dim_source>::Identity(), 
        const Vector<dim_target>& translation = Vector<dim_target>::Zero());
    
    // eg, for dim_source=2, this matrix acts on a vector of the form [x y 1]^t
    MatrixRect<1+dim_target, 1+dim_source> single_matrix_representation() const;

    // evaluation
    Vector<dim_target> operator()(const Vector<dim_source>& p) const;

    // composition (A*B)(x) := A(B(x))
    template<uint64_t dim_presource>
    Affinity<dim_target, dim_presource> operator*(const Affinity<dim_source, dim_presource>& other) const;

    MatrixRect<dim_target, dim_source> linear_part;
    Vector<dim_target> translation_part;
};


template<uint64_t dim>
struct AffinityIso : public Affinity<dim, dim>
{
    AffinityIso(const Matrix<dim>& linear = Matrix<dim>::Identity(), const Vector<dim>& translation = Vector<dim>::Zero());
    AffinityIso(const std::vector<Vector<dim>>& source, const std::vector<Vector<dim>>& target); // procrustes
    AffinityIso(const Affinity<dim, dim>& other);

    // this overload ensures that the composition of two AffinityIso objects is another AffinityIso object (not merely an Affinity object)
    AffinityIso<dim> operator*(const AffinityIso<dim>& other) const;

    AffinityIso<dim> inverse() const;

    // todo: implement this
    // return false if and only if no fixed points exist
    // the set of fixed points have the form particular_fixed_point + linear combo of basis elements
    bool fixed_points(Vector<dim>& particular_fixed_point, std::vector<Vector<dim>>& basis) const;
};


#include "affinity.hpp"