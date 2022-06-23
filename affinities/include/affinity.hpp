
#include "misc_math.h"

template<uint64_t dim_target, uint64_t dim_source>
Affinity<dim_target, dim_source>::Affinity(const MatrixRect<dim_target, dim_source>& linear, const Vector<dim_target>& translation) :
    linear_part(linear), translation_part(translation)
{
}


template<uint64_t dim_target, uint64_t dim_source>
MatrixRect<1+dim_target, 1+dim_source> Affinity<dim_target, dim_source>::single_matrix_representation() const
{
    MatrixRect<1+dim_target, 1+dim_source> rtn = MatrixRect<1+dim_target, 1+dim_source>::Zero();
    rtn.topLeftCorner(dim_target, dim_source) = linear_part;
    rtn.topRightCorner(dim_target, 1) = translation_part;
    rtn(dim_target, dim_source) = 1;

    return rtn;
}

template<uint64_t dim_target, uint64_t dim_source>
Vector<dim_target> Affinity<dim_target, dim_source>::operator()(const Vector<dim_source>& p) const
{
    return linear_part*p + translation_part;
}

template<uint64_t dim_target, uint64_t dim_source>
template<uint64_t dim_presource>
Affinity<dim_target, dim_presource> Affinity<dim_target, dim_source>::operator*(const Affinity<dim_target, dim_source>& other) const
{
    auto linear = linear_part*other.linear_part;
    auto trans = linear_part*other.translation_part + translation_part;
    return Affinity<dim_target, dim_presource>(linear, trans);
}



template<uint64_t dim>
AffinityIso<dim>::AffinityIso(const Matrix<dim>& linear, const Vector<dim>& translation) : 
    Affinity<dim, dim>(linear, translation)
{
    assert(linear.determinant() != 0);
}


template<uint64_t dim>
AffinityIso<dim> procrustes(const std::vector<Vector<dim>>& source, const std::vector<Vector<dim>>& target, bool force_special_orthogonal = false);


template<uint64_t dim>
AffinityIso<dim>::AffinityIso(const Affinity<dim, dim>& other) : AffinityIso(other.linear_part, other.translation_part)
{}

template<uint64_t dim>
AffinityIso<dim> AffinityIso<dim>::operator*(const AffinityIso<dim>& other) const
{
    return AffinityIso<dim>((*this)*other);
}


template<uint64_t dim>
AffinityIso<dim> AffinityIso<dim>::inverse() const
{
    auto linear = Affinity<dim, dim>::linear_part.inverse();
    auto translation = -linear*(Affinity<dim, dim>::translation_part);
    return AffinityIso<dim>(linear, translation);
}

