#include "conic.h"
#include "misc_math.h"
#include "misc.h"

Conic::Conic() : Conic(0,0,0,0,0,0)
{
    Matrix<3> H;
    get_coefficients(H);
    classification = Classification(H);
}
Conic::Conic(double a_, double b_, double c_, double d_, double e_, double f_) : a(a_), b(b_), c(c_), d(d_), e(e_), f(f_)
{
    Matrix<3> H;
    get_coefficients(H);
    classification = Classification(H);
}
Conic::Conic(const Vector<6>& q) : a(q(0)), b(q(1)), c(q(2)), d(q(3)), e(q(4)), f(q(5))
{
    Matrix<3> H;
    get_coefficients(H);
    classification = Classification(H);
}
Conic::Conic(const Matrix<3>& K) : a(K(0,0)), b(K(1,0) + K(0,1)), c(K(1,1)), d(K(2,0) + K(0,2)), e(K(2,1) + K(1,2)), f(K(2,2))
{
    Matrix<3> H;
    get_coefficients(H);
    classification = Classification(H);
}

void Conic::get_coefficients(Vector<6>& coefficients) const
{
    coefficients(0) = a;
    coefficients(1) = b;
    coefficients(2) = c;
    coefficients(3) = d;
    coefficients(4) = e;
    coefficients(5) = f;
}

void Conic::get_coefficients(Matrix<3>& H) const
{
    H(0,0) = a;
    H(1,0) = H(0,1) = b/2;
    H(1,1) = c;
    H(2,0) = H(0,2) = d/2;
    H(2,1) = H(1,2) = e/2;
    H(2,2) = f;
}

Conic Conic::transform(const Matrix<3>& B) const
{
    Matrix<3> A;
    get_coefficients(A);
    Matrix<3> temp = B.transpose()*A*B;
    return Conic(temp);
}

std::string Conic::to_string() const
{
    std::string rtn;
    rtn += std::to_string(a) + "*x^2 + ";
    rtn += std::to_string(b) + "*xy + ";
    rtn += std::to_string(c) + "*y^2 + ";
    rtn += std::to_string(d) + "*x + ";
    rtn += std::to_string(e) + "*y + ";
    rtn += std::to_string(f) + " = 0";
    return rtn;
}

Conic Conic::best_fit_conic(const std::vector<Vector<2>>& points) // performs a least squares fit for the coefficients
{
    Eigen::MatrixXd s = -1*Eigen::MatrixXd::Ones(points.size(),1);  

    Eigen::MatrixXd S(points.size(), 5);                            // sample data
    for (uint64_t i=0; i<points.size(); ++i)            
    {
        auto x = points[i](0);
        auto y = points[i](1);
        S(i,0) = x*x;
        S(i,1) = x*y;
        S(i,2) = y*y;
        S(i,3) = x;
        S(i,4) = y;
    }

    auto T = S.transpose()*S;
    auto t = S.transpose()*s;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(T);
    auto x = solver.solve(t);

    return Conic(x(0), x(1), x(2), x(3), x(4), 1);
}

void Conic::write_sampled_ellipse_to_vtk_file(const Conic& conic, const std::string& filename)
{
    const uint64_t num_points = 64; // todo: put this magic number somewhere else
    assert(num_points > 0);

    if (conic.classification.is_degenerate())
    {
        std::cerr << "currently, only handling nondegenerate conics in function write_sampled_ellipse_to_vtk_file" << std::endl;
        exit(1);
    }

    if (!conic.classification.is_ellipse())
    {
        std::cerr << "currently, only handling ellipses in function write_sampled_ellipse_to_vtk_file" << std::endl;
        exit(1);
    }

    Affinity<2,2> aligning_transformation = conic.axis_align();

    Conic aligned_conic = conic.transform(aligning_transformation.single_matrix_representation());

    if ((aligned_conic.classification).is_degenerate() || !aligned_conic.classification.is_ellipse())
    {
        std::cerr << "should not get here!!" << std::endl;
        exit(1);
    }

    if (!aligned_conic.isAxisAligned(1e-5))
    {
        std::cerr << "axis aligned failed" << std::endl;
        exit(1);
    }

    if (fabs(aligned_conic.f) < 1e-8 || fabs(aligned_conic.a) < 1e-8 || fabs(aligned_conic.c) < 1e-8)
    {
        std::cerr << std::endl << "conic: " << conic.to_string() << std::endl;
        std::cerr << std::endl << "aligned_conic: " << aligned_conic.to_string() << std::endl;
        std::cerr << "conic classification: " << conic.classification.to_string() << std::endl;
        std::cerr << "aligned conic classification: " << aligned_conic.classification.to_string() << std::endl;
        Matrix<3> coefficients;
        conic.get_coefficients(coefficients);
        std::cerr << "conic coefficients:" << coefficients << std::endl;
        aligned_conic.get_coefficients(coefficients);
        std::cerr << "aligned conic coefficients:" << coefficients << std::endl;

        exit(1);
    }

    // the aligned conic is axx + cyy + f = 0 where a and c have the same sign but not f
    if (aligned_conic.a * aligned_conic.c <= 0 || aligned_conic.a * aligned_conic.f >= 0)
    {
        std::cerr << "sign check failed" << std::endl;
        std::cerr << "aligned_conic: " << aligned_conic.to_string() << std::endl;
        exit(1);
    }

    // sample points along the axis aligned conic and then undo the aligning transformation

    Affinity<2,2> undo_aligning_transformation = aligning_transformation.inverse();

    std::vector<Vector<2>> samples;
    const double dtheta = 2*M_PI / num_points;
    const double alpha = sqrt(-aligned_conic.f/aligned_conic.a);
    const double beta = sqrt(-aligned_conic.f/aligned_conic.c);

    for (uint64_t i=0; i<num_points; ++i)
    {
        double theta = i*dtheta;
        Vector<2> sample;
        sample(0) = alpha*cos(theta);
        sample(1) = beta*sin(theta);
        samples.push_back(undo_aligning_transformation(sample));
    }

    Miscellaneous::write_points_to_vtk_file<2>(samples, filename, true);
}

Affinity<2,2> Conic::axis_align() const
{
    Matrix<3> H;
    get_coefficients(H);
    auto H0 = H.topLeftCorner(2,2);
    auto b = H.topRightCorner(2,1);

    // note: H0 is real symmetric (ie self adjoint) so it is diagonalizable and it diagonalizing matrix is orthogonal
    Eigen::SelfAdjointEigenSolver<Matrix<2>> es(H0);
    Matrix<2> V = es.eigenvectors();
    assert(V.isUnitary());

    // to get the translation part of the affinity, we approximate a solution to V*(H0 * x + b) = 0

    Eigen::ColPivHouseholderQR<Matrix<2>> solver(V*H0);
    Vector<2> x = solver.solve(-V*b);

    if (!(V*H0*x).isApprox(-V*b))
        std::cerr << "failed to find axis aligning transformation for the conic" << std::endl;

    return Affinity<2,2>(V, x);
}

bool Conic::isAxisAligned(const double tolerance) const
{
    Matrix<3> H;
    get_coefficients(H);
    for (uint64_t i=0; i<3; ++i)
        for (uint64_t j=0; j<3; ++j)
            if (i!=j && fabs(H(i,j))>tolerance)
                return false;

    return true;
}

std::array<double, 2> Conic::semimajor_and_semiminor_axes() const
{
    assert(!classification.is_degenerate());
    
    Affinity<2,2> aligning_transformation = axis_align();
    Conic conic0 = transform(aligning_transformation.single_matrix_representation());

    double a1 = sqrt(std::max<double>(1/fabs(conic0.a), 1/fabs(conic0.c)));
    double b1 = sqrt(std::min<double>(1/fabs(conic0.a), 1/fabs(conic0.c)));

    return std::array<double, 2>({a1,b1});
}

double Conic::eccentricity() const
{
    if (classification.is_degenerate())
    {
        std::cerr << "error:  cannot compute eccentricity of a degenerate conic"  << std::endl;
        assert(false);
    }

    if (classification.is_parabola())
    {
        return 1;
    }
    
    Matrix<3> H;
    get_coefficients(H);            

    double det = H.determinant();
    assert(det != 0); // since conic is nondegenerate
    double mu = det > 0 ? -1 : 1;

    double disc = (a-c)*(a-c) + b*b;
    double disc_sqrt = sqrt(disc);

    double num = 2*disc_sqrt;
    double denom = mu*(a+c) + disc_sqrt;
    assert(denom != 0);

    double ecc_sq = num/denom;
    assert(ecc_sq >= 0);

    double ecc = sqrt(ecc_sq);

    return ecc;
}

bool Conic::verify_eccentricity(double ecc, double tolerance) const
{
    double del = a*c-b*b/4;
    double temp = (a+c)*(a+c) - 4*del;
    double ecc2 = ecc*ecc;
    return fabs(del*ecc2*ecc2 + temp*(ecc2 - 1)) <= tolerance;
}

Conic::Classification::Classification(const Matrix<3>& Q)
{

    if (tolerance <= 0)
    {
        std::cerr << "invalide tolerance" << std::endl;
        exit(1);
    }

    double detQ = Q.determinant();

    Matrix<2> M = Q.topLeftCorner(2,2);
    double detM = M.determinant();

    if (fabs(detQ) > tolerance) // non-degenerate
    {
        if (detM > tolerance)
            code = ClassificationCode::NON_DEGENERATE_ELLIPSE;
        else if (detM < -tolerance)
            code = ClassificationCode::NON_DEGENERATE_HYPERBOLA;
        else
            code = ClassificationCode::NON_DEGENERATE_PARABOLA;
    }
    else
    {
        if (detQ > tolerance)
            code = ClassificationCode::SINGLE_POINT;
        else if (detQ < -tolerance)
            code = ClassificationCode::OBLIQUE_LINES;
        else
        {
            auto a = Q(0,0);
            auto c = Q(1,1);
            auto d = Q(0,1);
            auto e = Q(0,2);
            auto f = Q(2,2);

            const double l = d*d + e*e;
            const double r = (a+c)*f;
            if (l>r + tolerance)
                code = ClassificationCode::PARALLEL_LINES_DISTINCT_REAL;
            else if(l<r - tolerance)
                code = ClassificationCode::PARALLEL_LINES_IMAGINARY;
            else
                code = ClassificationCode::PARALLEL_LINES_COINCIDENT_AND_REAL;
        }
    }
}

std::string Conic::Classification::to_string() const
{
    switch (code)
    {
    case Conic::Classification::NON_DEGENERATE_ELLIPSE:
        return "NON_DEGENERATE_ELLIPSE";
    case Conic::Classification::NON_DEGENERATE_HYPERBOLA:
        return "NON_DEGENERATE_HYPERBOLA";
    case Conic::Classification::NON_DEGENERATE_PARABOLA:
        return "NON_DEGENERATE_PARABOLA";

    case Conic::Classification::OBLIQUE_LINES:
        return "OBLIQUE_LINES";

    case Conic::Classification::PARALLEL_LINES_DISTINCT_REAL:
        return "PARALLEL_LINES_DISTINCT_REAL";
    case Conic::Classification::PARALLEL_LINES_COINCIDENT_AND_REAL:
        return "PARALLEL_LINES_COINCIDENT_AND_REAL";
    case Conic::Classification::PARALLEL_LINES_IMAGINARY:
        return "PARALLEL_LINES_IMAGINARY";

    case Conic::Classification::SINGLE_POINT:
        return "SINGLE_POINT";
    }

    return "whaaaaaaat???";
}

bool Conic::Classification::is_degenerate() const
{
    return 
    code != ClassificationCode::NON_DEGENERATE_ELLIPSE &&
    code != ClassificationCode::NON_DEGENERATE_HYPERBOLA && 
    code != ClassificationCode::NON_DEGENERATE_PARABOLA;
}
bool Conic::Classification::is_parabola() const
{
    return 
    code == ClassificationCode::NON_DEGENERATE_PARABOLA ||
    code == ClassificationCode::PARALLEL_LINES_DISTINCT_REAL ||
    code == ClassificationCode::PARALLEL_LINES_COINCIDENT_AND_REAL ||
    code == ClassificationCode::PARALLEL_LINES_IMAGINARY;
}
bool Conic::Classification::is_hyperbola() const
{
    return 
    code == ClassificationCode::NON_DEGENERATE_HYPERBOLA ||
    code == ClassificationCode::OBLIQUE_LINES;
}
bool Conic::Classification::is_ellipse() const
{
    return 
    code == ClassificationCode::NON_DEGENERATE_ELLIPSE ||
    code == ClassificationCode::SINGLE_POINT;
}


