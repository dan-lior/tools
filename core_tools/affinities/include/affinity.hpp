
#include "misc_math.h"

template<uint64_t dim_target, uint64_t dim_source>
Affinity<dim_target, dim_source>::Affinity(const MatrixRect<dim_target, dim_source>& linear, const Vector<dim_target>& translation) :
    linear_part(linear), translation_part(translation)
{
}

template<uint64_t dim_target, uint64_t dim_source>
bool Affinity<dim_target, dim_source>::isInvertible() const
{
    return (dim_target == dim_source) && linear_part.determinant() != 0;
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
Affinity<dim_target, dim_presource> Affinity<dim_target, dim_source>::operator*(const Affinity<dim_source, dim_presource>& other) const
{
    auto linear = linear_part*other.linear_part;
    auto trans = linear_part*other.translation_part + translation_part;
    return Affinity<dim_target, dim_presource>(linear, trans);
}

template<uint64_t dim_target, uint64_t dim_source>
Affinity<dim_source, dim_target> Affinity<dim_target, dim_source>::inverse() const
{
    assert(isInvertible());

    auto l = linear_part.inverse();
    auto t = -l*translation_part;
    return Affinity<dim_source, dim_target>(l,t);
}


