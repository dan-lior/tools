#pragma once

#include "common_defs.h"
#include "misc_math.h"
#include "misc.h"

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
    Affinity<dim_target, dim_presource> operator*(const Affinity<dim_source, dim_presource>& other) const;

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

    static Affinity<dim_target, dim_source> procrustes(const std::vector<Vector<dim_source>>& source, const std::vector<Vector<dim_target>>& target)
    {
        static_assert(dim_target = dim_source, "procrustes only implemented for equal template parameters");
        constexpr uint64_t dim = dim_target;

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

    //    Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);

        Eigen::MatrixXd U = svd.matrixU();
        Eigen::MatrixXd V = svd.matrixV();

        Matrix<dim> R = V*U.transpose();

        // to force special orthogonal (ie det =1 and not det = -1)
        // if (R.determinant() < 0)
        // {
        //     std::cerr << "R.determinant(): " << R.determinant()<< std::endl;
        //     Matrix<dim> flip = Matrix<dim>::Identity();
        //     flip(dim-1, dim-1) = -1; 
        //     R = V * flip * U.transpose();
        //     std::cerr << "R.determinant(): " << R.determinant()<< std::endl;
        // }

        Affinity<dim, dim> AffinityIso1(Matrix<dim>::Identity(), -centroid_source);
        Affinity<dim, dim> AffinityIso2(R, Vector<dim>::Zero());
        Affinity<dim, dim> AffinityIso3(Matrix<dim>::Identity(), centroid_target);

        return AffinityIso3 * AffinityIso2 * AffinityIso1;
    }




    MatrixRect<dim_target, dim_source> linear_part;
    Vector<dim_target> translation_part;
};

#include "affinity.hpp"