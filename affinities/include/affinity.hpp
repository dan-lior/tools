
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
Affinity<dim_target, dim_presource> operator*(const Affinity<dim_source, dim_presource>& other) const
{
    return Affinity<dim_presource, dim_source>(linear_part*other.linear_part, linear_part*other.translation_part + translation_part);
}

template<uint64_t dim>
AffinityIso<dim>::AffinityIso(const Matrix<dim>& linear, const Vector<dim>& translation) : 
    Affinity<dim, dim>(linear, translation)
{
    assert(linear.determinant() != 0);
}

template<uint64_t dim>
AffinityIso<dim>::AffinityIso(const std::vector<Vector<dim>>& source, const std::vector<Vector<dim>>& target, bool force_special_orthogonal = false) : 
    Affinity<dim, dim>(source, target)
{

    assert(!source.empty());
    assert(source.size() == target.size());
    const uint64_t n = source.size();

    const Vector<dim> centroid_source = MiscellaneousMath::centroid<dim>(source);
    const Vector<dim> centroid_target = MiscellaneousMath::centroid<dim>(target);
        

    Eigen::MatrixXd A(dim, n);
    Eigen::MatrixXd B(dim, n);
    for (uint64_t i=0; i<n; ++i)
    {
        A.col(i) = source[i] - centroid_source;
        B.col(i) = target[i] - centroid_target;
    }

    Eigen::MatrixXd M = A*B.transpose();

//    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M); // this is the newer version specified in the docs, but it doesn't seem to build
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV); // this depracated version still works

    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();

    Matrix<dim> R = V*U.transpose();

    if (force_special_orthogonal && R.determinant() < 0)
    {
        Matrix<dim> flip = Matrix<dim>::Identity();
        flip(dim-1, dim-1) = -1; 
        R = V * flip * U.transpose();
    }

    AffinityIso<dim> AffinityIso1(Matrix<dim>::Identity(), -centroid_source);
    AffinityIso<dim> AffinityIso2(R, Vector<dim>::Zero());
    AffinityIso<dim> AffinityIso3(Matrix<dim>::Identity(), centroid_target);

    return AffinityIso3 * AffinityIso2 * AffinityIso1;
}

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

