#include "polyline.h"
#include "misc.h"

namespace
{
    double registration_error(std::vector<Vector<3>>& as, std::vector<Vector<3>>& bs)
    {
        double sum_squared_error = 0;
        assert(as.size() == bs.size());

        for (uint64_t i=0; i< as.size(); ++i)
            sum_squared_error += (as[i] - bs[i]).squaredNorm();

        return sqrt(sum_squared_error / as.size());
    }

    // returns a transformation that best aligns floating with fixed
    Affinity<3,3> single_icp_iterate(const Polyline& fixed, const Polyline& floating, bool vanilla)
    {
        std::vector<Vector<3>> targets;

        if (vanilla)
            for (auto q : floating.vertices)
            {
                uint64_t i = MiscellaneousMath::closest_point<3>(fixed.vertices, q);
                targets.push_back(fixed.vertices[i]);
            }
        else
            for (auto q : floating.vertices)
            {
                uint64_t i = MiscellaneousMath::closest_point<3>(fixed.vertices, q);
                targets.push_back(MiscellaneousMath::projection_to_line<3>(fixed.vertices[i], fixed.tangent_directions[i], q));
            }

        // todo: fix this
        //return Affinity<3,3>::procrustes(floating.vertices, targets);
        return Affinity<3,3>();
    }

    // returns the angle between (p2-p1) and (p3-p2)
    double external_angle(const Vector<3>& p1, const Vector<3>& p2, const Vector<3>& p3)
    {
        const Vector<3> d = p2-p1;
        const Vector<3> e = p3-p2;
        assert(d.norm() * e.norm() > 0);

        const double cos_phi = (p2-p1).dot(p3-p2) / (d.norm() * e.norm());
        const double phi = std::acos(cos_phi);

        return phi;
    }

    // returns approximation to curvature at p2
    double three_point_curvature(const Vector<3>& p1, const Vector<3>& p2, const Vector<3>& p3)
    {
        // from Max Planck Institute research report:
        // "Analysis and Design of Discrete Normals and Curvatures"
        // Langer, Belyaev, Seidel
        // March 2005

        const double phi = external_angle(p1,p2,p3);

        double d = (p2-p1).norm();
        double e = (p3-p2).norm();
        assert(d*e > 0);

        const double k = 2*phi / (d + e);

        return k;
    }

    // returns approximation to curvature at p3
    double five_point_curvature(const Vector<3>& p1, const Vector<3>& p2, const Vector<3>& p3, 
        const Vector<3>& p4, const Vector<3>& p5)
    {
        // from Max Planck Institute research report:
        // "Analysis and Design of Discrete Normals and Curvatures"
        // Langer, Belyaev, Seidel
        // March 2005

        const double phi_2 = external_angle(p1, p2, p3);
        const double phi_4 = external_angle(p3, p4, p5);

        const double c = (p2-p1).norm();
        const double d = (p3-p2).norm();
        const double e = (p4-p3).norm();
        const double f = (p5-p4).norm();

        const double k2 = 2*phi_2 / (c+d);
        const double k4 = 2*phi_4 / (e+f);

        const double k = ((2*e+f)*k2 + (c+2*d)*k4) / (c + 2*(d+e) + f);

        return k;
    }

    // returns approximation to torsion at p3
    double five_point_torsion(const Vector<3>& p1, const Vector<3>& p2, const Vector<3>& p3, 
        const Vector<3>& p4, const Vector<3>& p5, const Vector<3>& t3)
    {
        // from Max Planck Institute research report:
        // "Analysis and Design of Discrete Normals and Curvatures"
        // Langer, Belyaev, Seidel
        // March 2005

        const double k2 = three_point_curvature(p1, p2, p3);
        const double k4 = three_point_curvature(p3, p4, p5);

        const Vector<3> c_ = (p2-p1);
        const Vector<3> d_ = (p3-p2);
        const Vector<3> e_ = (p4-p3);
        const Vector<3> f_ = (p5-p4);

        const double c = c_.norm();
        const double d = d_.norm();
        const double e = e_.norm();
        const double f = f_.norm();

        const double denom1 = c + 2*(d+e) + f;
        assert(denom1 > 0);

        const double k = ((2*e+f)*k2 + (c+2*d)*k4) / denom1;
        assert(k != 0);
        const double k_prime = 3 * (k4  - k2) / denom1;

        const Vector<3> b2 = (c_.cross(d_)).normalized();
        const Vector<3> b3 = (d_.cross(e_)).normalized();
        const Vector<3> b4 = (e_.cross(f_)).normalized();

        const double mu_d = (b4.cross(b3)).dot(t3);
        const double mu_e = (b3.cross(b2)).dot(t3); // possibly negate this

        const double tau_d = 3*mu_d / (c+d+e);
        const double tau_e = 3*mu_e / (d+e+f);
        
        const double num2 = (f+2*e-d)*tau_d + (c+2*d-e)*tau_e;
        const double denom2 = c+d+e+f;
        const double tau_hat = num2/denom2;

        const double num3 = c*e - d*f +d*d - e*e;
        const double denom3 = 3*denom2;
        const double tau = tau_hat * (1 + (k_prime / k) * num3 / denom3);

        return tau;
    }
}

    Polyline::Polyline(const std::vector<Vector<3>>& vertices_, uint64_t smoothing_radius_) 
        : vertices(vertices_), smoothing_radius(smoothing_radius_)
    {   
        assert(vertices.size() > 1);
        assert(smoothing_radius > 0);
        generate_smooth_tangents();
    }

    Polyline::Polyline(const std::string& filename, uint64_t smoothing_radius_): 
        Polyline(Miscellaneous::read_points_from_vtk_file<3>(filename), smoothing_radius_)
    {
    }

    void Polyline::generate_smooth_tangents() 
    {
        const uint64_t n = vertices.size();
        tangent_directions.resize(n);
        
        if (n == 2)
        {
            Vector<3> temp = (vertices[1]-vertices[0]).normalized();
            tangent_directions[0] = tangent_directions[1] = temp;
        }
        else if (n > 2)
        {
            assert(smoothing_radius >= 1);
            if (smoothing_radius == 1)
            {
                // first, compute the directions of the internal vertices (ie, all but the endpoints)
                for (uint64_t i = 1; i<n-1; ++i)
                {
                    // from Max Planck Institute research report:
                    // "Analysis and Design of Discrete Normals and Curvatures"
                    // Langer, Belyaev, Seidel
                    // March 2005

                    const Vector<3> d_ = (vertices[i] - vertices[i-1]);
                    const Vector<3> e_ = (vertices[i+1] - vertices[i]);
                    const double d = d_.norm();
                    const double e = e_.norm();

                    tangent_directions[i] = d*e/(d+e)*(d_/(d*d) + e_/(e*e));;
                }

                // second, assign the directions of the endpoints by constant extension

                tangent_directions[0] = tangent_directions[1];
                tangent_directions[n-1] = tangent_directions[n-2];
            }
            else
            {
                assert(smoothing_radius >= 2);

                // for each vertex v of the polyline:
                //  - select an iterval of vertices of given radius and centered on v
                //  - fit a line to these samples and store this "tangent direction at v"
                        
                for (uint64_t i = 0; i<vertices.size(); ++i)
                {
                    std::vector<Vector<3>> points;
                    const uint64_t i_begin = i < smoothing_radius ? 0 : i-smoothing_radius;
                    const uint64_t i_end = std::min<uint64_t>(vertices.size(), i+smoothing_radius+1);
                    for (uint64_t j = i_begin; j<i_end; ++j)
                        points.push_back(vertices[j]);

                    auto line = MiscellaneousMath::fit_line3d(points);
                    tangent_directions[i] = line[1];
                }
            }
        }
    }

    std::vector<double> Polyline::distances_from_start() const
    {
        std::vector<double> distances;
        distances.push_back(0);
        for (uint64_t i=1; i<vertices.size(); ++i)
        {
            double interval_length = (vertices[i] - vertices[i-1]).norm();
            distances.push_back(interval_length + distances.back());
        }
        return distances;
    }

    Polyline Polyline::reverse_direction() const
    {
        std::vector<Vector<3>> vertices_reversed = vertices;
        std::reverse(vertices_reversed.begin(), vertices_reversed.end());
        return Polyline(vertices_reversed, smoothing_radius);
    }

    void Polyline::display_stats() const
    {
        std::vector<double> interval_lengths;
        double sum=0;
        for (uint64_t i=1; i<vertices.size(); ++i)
        {
            double interval_length = (vertices[i] - vertices[i-1]).norm();
            interval_lengths.push_back(interval_length);
            sum +=  interval_length;
        }
        std::cerr << "number of intervals: " << interval_lengths.size() << std::endl;
        std::cerr << "first interval length: " << interval_lengths.front() << std::endl;
        std::cerr << "last interval length: " << interval_lengths.back() << std::endl;
        
        std::cerr << "min interval length: " << *std::min_element(interval_lengths.begin(), interval_lengths.end()) << std::endl;
        std::cerr << "max interval length: " << *std::max_element(interval_lengths.begin(), interval_lengths.end()) << std::endl;
        std::cerr << "avg interval length: " << sum / interval_lengths.size() << std::endl;
    }

    Polyline Polyline::resample(double mm) const
    {
        assert(mm > 0);
        const double length = distances_from_start().back();
        const uint64_t n = 1 + std::ceil(length / mm); // number of uniform vertices
        const double d = length / (n-1);
        // d approximates mm and divides the length of the polyline evenly

        auto distance_from_start = distances_from_start();
        assert(distance_from_start.size() == vertices.size());

        std::vector<Vector<3>> vertices_uniform;
        vertices_uniform.push_back(vertices[0]);

        for (uint64_t i=0; i<vertices.size(); ++i)
            while (d*vertices_uniform.size() <= distance_from_start[i])
            {
                auto temp = vertices[i] - vertices_uniform.back();
                auto u = vertices_uniform.back() + (d/temp.norm())*temp;
                vertices_uniform.push_back(u);
            }

        // assert(vertices_uniform.size() == n);
        int disc = vertices_uniform.size() - n;
        assert(disc <=1 && disc >=-1);

        return Polyline(vertices_uniform, smoothing_radius);
    }

    Polyline Polyline::reset_smoothing_radius(uint64_t new_smoothing_radius) const
    {
        return Polyline(vertices, new_smoothing_radius);
    }
    
    void Polyline::write_to_file(const std::string& filename) const
    {
        Miscellaneous::write_points_to_vtk_file<3>(vertices, filename);
    }

    std::vector<double> Polyline::curvature() const
    {
        const uint64_t n = vertices.size();
        std::vector<double> curvatures(n, 0);

        if (n >= 3)
        {
            // convention for first two and last two vertices
            curvatures[0] = curvatures[1] = three_point_curvature(vertices[0], vertices[1], vertices[2]);
            curvatures[n-1] = curvatures[n-2] = three_point_curvature(vertices[n-3], vertices[n-2], vertices[n-1]);

            // compute curvatures at the internal vertices
            for (uint64_t i=2; i<n-2; ++i)
                curvatures[i] = five_point_curvature(vertices[i-2], vertices[i-1], vertices[i], vertices[i+1], vertices[i+2]);
        }

        return curvatures;
    }

    std::vector<double> Polyline::torsion() const
    {
        const uint64_t n = vertices.size();
        std::vector<double> torsion(n, 0);

        if (n >= 5)
        {
            // compute torsion at the internal vertices
            for (uint64_t i=3; i<n-3; ++i)
                torsion[i] = five_point_torsion(vertices[i-2], vertices[i-1], vertices[i], vertices[i+1], vertices[i+2], tangent_directions[i]);

            // convention for first two and last two vertices
            torsion[0] = torsion[1] = torsion[2];
            torsion[n-1] = torsion[n-2] = torsion[n-3];
        }

        return torsion;
    }

    Polyline Polyline::transform(const Affinity<3,3>& affinity) const
    {
        std::vector<Vector<3>> vertices_transformed;
        for (auto v : vertices)
            vertices_transformed.push_back(affinity(v));

        return Polyline(vertices_transformed, smoothing_radius);
    }

    Polyline Polyline::gaussian_perturbation(double mu, double sigma) const
    {
        std::vector<Vector<3>> vertices_transformed;
        for (auto v : vertices)
            vertices_transformed.push_back(v + MiscellaneousMath::gaussian_vector<3>(mu, sigma));

        return Polyline(vertices_transformed, smoothing_radius);

    }

    // Affinity<3,3> Polyline::register_icp(const Polyline& fixed, const Polyline& floating, double convergence_tolerance, uint64_t max_iterations, bool vanilla)
    // {
    //     Affinity<3,3> A; // identity
    //     for (uint64_t i = 0; i < max_iterations; ++i)
    //     {
    //         Affinity<3,3> A0 = single_icp_iterate(fixed, floating.transform(A), vanilla);            
    //         const double norm_A0 = (A0.single_matrix_representation() - Matrix<4>::Identity()).norm();
    //         std::cerr << "norm A0: " << norm_A0 << std::endl;
            
    //         A = A0*A;
    //         if (norm_A0 < convergence_tolerance)
    //         {
    //             std::cerr << "icp converged" << std::endl;
    //             return A;
    //         }
    //     }

    //     std::cerr << "icp failed to converge" << std::endl;
    //     return A;
    // }
