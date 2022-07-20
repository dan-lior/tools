#pragma once

#include "common_defs.h"

namespace MiscellaneousMath
{

    template <uint64_t dim>
    Vector<dim> gaussian_vector(double mu, double sigma);

    template<uint64_t dim>
    uint64_t closest_point(const std::vector<Vector<dim>>& points, const Vector<dim>& query_point);

    template<uint64_t dim>
    Vector<dim> projection_to_line(const Vector<dim>& point_on_line, const Vector<dim>& direction_of_line, const Vector<dim>& query_point);

    template<uint64_t dim>
    std::array<Vector<dim>, 2> bounding_box(const std::vector<Vector<dim>>& points);

    template<typename T>
    std::vector<std::vector<T>> cartesian_product(const std::vector<std::vector<T>>& sequences);

    // TODO: implement this directly (will result in simpler, faster code)
    template<uint64_t dim>
    std::vector<Index<dim>> cartesian_product(const Index<dim>& min_corner, const Index<dim>& max_corner);

    std::array<Vector<3>,2> fit_line3d(const std::vector<Vector<3>>& points); 

    // returns a matrix that projects to the hyperplane perpendicular to n
    template<uint64_t dim>
    MatrixRect<dim-1, dim> projection_matrix(const Vector<dim>& normal);

    Matrix<3> rotation(const Vector<3>& axis, double angle_in_radians);

    double sample_variance(const std::vector<double>& xs);

    double mean(const std::vector<double>& xs);

    template<uint64_t dim>
    std::pair<uint64_t, std::array<uint64_t,3>> rank_and_signature(Matrix<dim> A, double tolerance);

    template<uint64_t dim>
    double diameter(const std::vector<Vector<dim>>& points);

    template<uint64_t dim>
    Vector<dim> centroid(const std::vector<Vector<dim>>& points);

    // returns the distinct real roots of x^2+bx+c = 0
    std::vector<double> real_quadratic_roots(double b, double c);

    // no "cardinal point" lies strictly between any two points
    // cardinal points are definitely on the convex hull of the set of points
    // if there are more than two distinct cardinal points, returns a positively oriented convex polygon
    // if there is exactly one cardinal point, the elements of points are identical
    // returns an empty container if points is empty
    std::vector<Vector<2>> distinct_cardinal_points(const std::vector<Vector<2>>& points);

    // if points are not collinear, returned points are orded in counter clockwise order and no three points are collinear
    std::vector<Vector<2>> convex_hull(const std::vector<Vector<2>>& points);

    // strict containment
    bool point_inside_convex_polygon(const Vector<2> p, const std::vector<Vector<2>>& convex_polygon);

    // returns negative value is the triangle is oriented clockwise
    double signed_area_of_triangle(const std::array<Vector<2>,3>& triangle);

    // returns true if size at least 3 and every 3 consecutive points form a positively oriented triangle
    bool is_convex_positively_oriented_polygon(const std::vector<Vector<2>>& polygon);

    // only certain that this gives the right answer for convex polygons ... look into this
    double area_of_polygon(const std::vector<Vector<2>>& polygon);

    std::array<double, 3> min_max_variance(const std::vector<Vector<2>>& points);

    std::array<double, 2> tilted_bbox_dims(const std::vector<Vector<2>>& points);


}

#include "misc_math.hpp"