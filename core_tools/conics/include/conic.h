    #pragma once

    #include "common_defs.h"
    #include "affinity.h"

    class Conic
    {
    public:        
        
        static void write_sampled_ellipse_to_vtk_file(const Conic& conic, const std::string& filename);

        Conic();
        Conic(double a, double b, double c, double d, double e, double f); // axx + bxy + cyy + dx + ey + f == 0
        Conic(const Vector<6>& q);
        Conic(const Matrix<3>& H);
        
        void get_coefficients(Vector<6>& coefficients) const;
        void get_coefficients(Matrix<3>& H) const;

        Conic transform(const Matrix<3>& B) const;

        Affinity<2,2> axis_align() const;

        bool isAxisAligned(const double tolerance) const;

        std::string to_string() const;

        std::array<double, 2> semimajor_and_semiminor_axes() const;

        double eccentricity() const;

        // performs a least squares fit for the coefficients ("algebraic fit")
        static Conic best_fit_conic(const std::vector<Vector<2>>& points); 

        bool verify_eccentricity(double ecc, double tolerance) const;

        struct Classification
        {
            static constexpr double tolerance = 1e-8;

            Classification(const Matrix<3>& Q = Matrix<3>::Zero());

            std::string to_string() const;

            bool is_degenerate() const;
            bool is_parabola() const;
            bool is_hyperbola() const;
            bool is_ellipse() const;

            enum ClassificationCode 
            {
                NON_DEGENERATE_ELLIPSE=0, NON_DEGENERATE_PARABOLA, NON_DEGENERATE_HYPERBOLA, 
                OBLIQUE_LINES,  // degenerate hyperbola
                PARALLEL_LINES_DISTINCT_REAL, PARALLEL_LINES_COINCIDENT_AND_REAL, PARALLEL_LINES_IMAGINARY, // degenerate parabola
                SINGLE_POINT // degenerate ellipse
            } code;

        } classification;

    private:

        double a;
        double b;
        double c;
        double d;
        double e;
        double f;
    };
