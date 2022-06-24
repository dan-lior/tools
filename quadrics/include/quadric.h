#pragma once

#include "affinity.h"
#include "common_defs.h"
#include "misc.h"
#include "misc_math.h"
#include "conic.h"

struct Quadric
{
    // coefficients represent the quadric:
    // a + bx + cy + dz + exx + fxy + gxz + hyy + iyz + jzz == 0

    Quadric();
    Quadric(double a_, double b_, double c_, double d_, double e_, double f_, double g_, double h_, double i_, double j_);
    Quadric(const Vector<10>& q);
    Quadric(const Matrix<4>& H);
        
    void get_coefficients(Vector<10>& coefficients) const;
    void get_coefficients(Matrix<4>& H) const
    {
        H(0,0) = e;
        H(1,1) = h;
        H(2,2) = j;
        H(3,3) = a;

        H(1,0) = H(0,1) = f/2;
        H(2,0) = H(0,2) = g/2;
        H(3,0) = H(0,3) = b/2;

        H(2,1) = H(1,2) = i/2;
        H(3,1) = H(1,3) = c/2;

        H(3,2) = H(2,3) = d/2;
    }

    // todo
    // Quadric transform(const Affinity<3>& affinity) const
    // {
    //     return Quadric();
    // }

    // Quadric transform(const Matrix<4>& projectivity) const
    // {
    //     return Quadric();
    // }

    std::string to_string() const
    {
        std::string rtn;
        rtn += std::to_string(a) + " + ";
        rtn += std::to_string(b) + "*x + ";
        rtn += std::to_string(c) + "*y + ";
        rtn += std::to_string(d) + "*z + ";
        rtn += std::to_string(e) + "*x^2 + ";
        rtn += std::to_string(f) + "*xy + ";
        rtn += std::to_string(g) + "*xz + ";
        rtn += std::to_string(h) + "*y^2 + ";
        rtn += std::to_string(i) + "*yz + ";
        rtn += std::to_string(j) + "*zz = 0";
        return rtn;
    }


    enum QuadricFitTypes { TYPE_GENERAL_QUADRIC = 0, TYPE_ROTSYM_QUADRIC, TYPE_PLANE, TYPE_SPHERE,                  // 0-3
        TYPE_GEN_CYL, TYPE_CIRC_CYL, TYPE_CONE, TYPE_CIRC_CONE,                                     // 4-7
        TYPE_ELLIPSOID_BIASED, TYPE_HYPERBOLOID_BIASED, TYPE_ELLIPSOID_OPT, TYPE_HYPERBOLOID_OPT,   // 8-11
        TYPE_HYPERBOLOID_1SHEET, TYPE_HYPERBOLOID_2SHEET, TYPE_PARABOLOID,                          // 12-14
        TYPE_PARABOLOID_ELLIPTICAL, TYPE_PARABOLOID_HYPERBOLIC,                                     // 15-16
        TYPE_ELL_CYL, TYPE_HYPER_CYL, TYPE_PARA_CYL,                                                // 17-19
        NUM_QUADRIC_TYPES } ;


    bool is_elliptic_cylinder(const double tolerance = 0.0)
    {
        std::vector<double> q;
        q.push_back(a);
        q.push_back(b);
        q.push_back(c);
        q.push_back(d);
        q.push_back(e);
        q.push_back(f);
        q.push_back(g);
        q.push_back(h);
        q.push_back(i);
        q.push_back(j);

        Matrix<3> Q;
        Q(0,0) = q[4];      //A
        Q(0,1) = q[5]/2;    //D
        Q(0,2) = q[6]/2;    //E 
        Q(1,0) = Q(0,1);
        Q(1,1) = q[7];      //B
        Q(1,2) = q[8]/2;    //F
        Q(2,0) = Q(0,2);
        Q(2,1) = Q(1,2);
        Q(2,2) = q[9];      //C

        Matrix<4> S;
        S.topLeftCorner<3,3>() = Q;
        S(0,3) = q[1]/2;    //G
        S(1,3) = q[2]/2;    //H
        S(2,3) = q[3]/2;    //I
        S(3,3) = q[0];      //J
        S(3,0) = S(0,3);
        S(3,1) = S(1,3);
        S(3,2) = S(2,3);

        auto rank_sig_Q = MiscellaneousMath::rank_and_signature<3>(Q, tolerance);
        auto rank_sig_S = MiscellaneousMath::rank_and_signature<4>(S, tolerance);

        if (rank_sig_Q.first != 2)
        {
            std::cerr << "Q has rank " << rank_sig_Q.first << " but 2 was expected" << std::endl;
            std::cerr << "S has rank " << rank_sig_S.first << std::endl;
            return false;
        } 
        if (rank_sig_S.first != 3)
        {
            std::cerr << "S has rank " << rank_sig_S.first << " but 3 was expected" << std::endl;
            std::cerr << "Q has rank " << rank_sig_Q.first << std::endl;
            return false;
        } 

        auto sig = rank_sig_Q.second;

        if (sig[0] == 2 && sig[1] == 0 && sig[2] == 1) return true;
        if (sig[0] == 0 && sig[1] == 2 && sig[2] == 1) return true;

        std::cerr << "Q has signature: " << sig[0] << " " << sig[1]  << " " << sig[2] << " not as expected" << std::endl;

        return false;
    }

    double eccentricity(double tolerance)
    {
        Matrix<4> H;
        get_coefficients(H);            

        // for a nondegenerate cylinder, exactly one of the eigenvalues of H is zero or near zero
        // the eigen vector associated to this zero eigenvalue is the axial direction of the cylinder
        // we project the surface along this axis to get a conic (an ellipse) and report its eccentricity.


        // convert to canonical form :
        Eigen::SelfAdjointEigenSolver<Matrix<4>> eigensolver(H);
        if (eigensolver.info() != Eigen::Success)
        {
            std::cerr << "eigensolver failed" << std::endl;
            assert(false);
        }

        // std::cerr << "matrix of eigenvectors: " << eigensolver.eigenvectors() << std::endl;


        double min_evalue = std::numeric_limits<double>::max();
        int imin_evalue=0;
        for (int i=0; i<4; ++i)
        {
            double ev = fabs(eigensolver.eigenvalues()[i]);
            if (ev < min_evalue)
            {
                min_evalue = ev;
                imin_evalue = i;
            }

        }
        
        assert(min_evalue < tolerance);

        auto temp = eigensolver.eigenvectors().col(imin_evalue);
        Vector<4> axis(temp(0), temp(1), temp(2), temp(3));

        // std::cerr << "imin_evalue: " << imin_evalue << std::endl;
        // std::cerr << "min_evalue: " << min_evalue << std::endl;
        // std::cerr << "temp: " << temp << std::endl;
        // std::cerr << "axis: " << axis << std::endl;



        auto Q = MiscellaneousMath::projection_matrix<4>(axis);
        // std::cerr << "projection matrix: " << std::endl;
        // std::cerr << Q << std::endl;

        Matrix<3> H0 = Q*H*Q.transpose();
        Conic conic(H0);

        assert(!conic.classification.is_degenerate());

        return conic.eccentricity();
    }

    static std::vector<Quadric> fit_all_types_of_quadrics(const std::string& filename);

    private:

    double a,b,c,d,e,f,g,h,i,j;
};

