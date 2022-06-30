#include "misc_math.h"

namespace
{
    bool top_three_make_left_turn(std::stack<Vector<2>> stk)
    {
        assert(stk.size() >= 3);
        std::array<Vector<2>,3> triangle;

        triangle[2] = stk.top();
        stk.pop();
        triangle[1] = stk.top();
        stk.pop();
        triangle[0] = stk.top();
        stk.push(triangle[1]);
        stk.push(triangle[2]);

        return MiscellaneousMath::signed_area_of_triangle(triangle) > 0;
    }

    void delete_one_below_top(std::stack<Vector<2>>& stk)
    {
        const uint64_t starting_size = stk.size();

        if (starting_size < 2 )
        {
            std::cerr << "logic error: in misc math" << std::endl;
            exit(1);
        }

        auto top = stk.top();
        stk.pop();

        if (stk.size() != starting_size -1 )
        {
            std::cerr << "logic error: in misc math" << std::endl;
            exit(1);
        }

        stk.pop();

        if (stk.size() != starting_size -2 )
        {
            std::cerr << "logic error: in misc math" << std::endl;
            exit(1);
        }

        stk.push(top);

        if (stk.size() != starting_size -1 )
        {
            std::cerr << "logic error: in misc math" << std::endl;
            exit(1);
        }

    }
}

Matrix<3> MiscellaneousMath::rotation(const Vector<3>& axis, double angle_in_radians)
{
    assert(axis.norm() > 0);

    // find an o.n. basis of perp space
    Eigen::JacobiSVD<Vector<3>> svd(axis, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto U = svd.matrixU();
    assert(U.isUnitary());
    // U rotates the x-axis to align with the given axis
        
    Matrix<3> R = Matrix<3>::Identity();
    double c = std::cos(angle_in_radians);
    double s = std::sin(angle_in_radians);
    R(1,1) = c;
    R(2,1) = s;
    R(1,2) = -s;
    R(2,2) = c;
    assert(R.isUnitary());
    assert(R.determinant() > 0);
    // R is a positive rotation by the given angle about the x-axis


    Matrix<3> rtn = U*R*(U.transpose());
    assert(rtn.isUnitary());
    return rtn;
}


std::array<Vector<3>,2> MiscellaneousMath::fit_line3d(const std::vector<Vector<3>>& points) 
{

    //verify that points contains at least two distinct elements
    assert(points.size() > 1);
    {
        bool at_least_two_distinct = false;
        for (auto p : points)
            if (p != points[0])
            {
                at_least_two_distinct = true;
                break;
            }
        if (!at_least_two_distinct)
        {
            std::cerr << "oops, not enough points to fit a line" << std::endl;
            exit(1);
        }
    }

    const Vector<3> centroid = MiscellaneousMath::centroid<3>(points);

    // power iteration algorithm (instead of svd)
    const double cos_angle_tolerance = 1.0e-6; // how to choose this appropriately? 
    // assert(cos_angle_tolerance > 0);
    // assert(cos_angle_tolerance < 1);
    double cos_angle = 0;

    Vector<3> direction = Vector<3>(0,-1,0); // any unit vector will do here
    while (fabs(1.0-cos_angle) > cos_angle_tolerance)
    {
        Vector<3> old_direction = direction;
        direction = Vector<3>::Zero();
        for (auto p : points)
        {
            const Vector<3> centered_point = p - centroid;
            direction = direction +  (centered_point.dot(old_direction) * centered_point);
        }

        // direction.normalize();
        direction = direction.normalized();
        
        cos_angle = direction.dot(old_direction);
    }

    Vector<3> rough_estimate = (points.back() - points.front()).normalized();
    if (rough_estimate.dot(direction) < 0)
        direction = -direction;


    std::array<Vector<3>,2> rtn;
    rtn[0] = centroid;
    rtn[1] = direction;
    return rtn;
}

double MiscellaneousMath::mean(const std::vector<double>& xs)
{
    assert(xs.size() > 1);

    double sum = 0;
    for (auto x:xs)
        sum += x;
        
    return sum / xs.size();
}

double MiscellaneousMath::sample_variance(const std::vector<double>& xs)
{
    const double x_mean = mean(xs);
    double sample_variance=0;

    for (auto x:xs)
    {
        double deviation = x-x_mean;
        sample_variance += deviation*deviation;
    }

    assert(xs.size() > 1);
    sample_variance /= (xs.size() - 1);

    return sample_variance;
}

std::vector<double> MiscellaneousMath::real_quadratic_roots(double b, double c)
{
    std::vector<double> roots;

    double disc = b*b - 4*c;

    if (disc == 0)
        roots.push_back(-b/2);
    else if (disc > 0)
    {
        roots.push_back((b>0 ? -b - sqrt(disc) : -b - sqrt(disc)) / 2);
        roots.push_back(c / roots.front());
    }

    return roots;
}

std::vector<Vector<2>> MiscellaneousMath::distinct_cardinal_points(const std::vector<Vector<2>>& points)
{
    if (points.empty())
        return std::vector<Vector<2>>();

    double minx = std::numeric_limits<double>::max();
    double miny = std::numeric_limits<double>::max();
    double maxx = -std::numeric_limits<double>::max();
    double maxy = -std::numeric_limits<double>::max();

    uint64_t i_minx = std::numeric_limits<uint64_t>::max(); // index of point with least x and then least y
    uint64_t i_miny = std::numeric_limits<uint64_t>::max(); // index of point with least y and then greatest x
    uint64_t i_maxx = std::numeric_limits<uint64_t>::max(); // index of point with greatest x and then greatest y
    uint64_t i_maxy = std::numeric_limits<uint64_t>::max(); // index of point with greatest y and then least x

    for (uint64_t i =0; i<points.size(); ++i)
    {
        const Vector<2>& p = points[i];

        if (p(0) < minx || (p(0) == minx && p(1) < points[i_minx](1)))
        {
            i_minx = i;
            minx = p(0);
        }

        if (p(1) < miny || (p(1) == miny && p(0) > points[i_miny](0)))
        {
            i_miny = i;
            miny = p(1);
        }

        if (p(0) > maxx || (p(0) == maxx && p(1) > points[i_maxx](1)))
        {
            i_maxx = i;
            maxx = p(0);
        }

        if (p(1) > maxy || (p(1) == maxy && p(0) < points[i_maxy](0)))
        {
            i_maxy = i;
            maxy = p(1);
        }
    }

    assert(i_minx < points.size());
    assert(i_miny < points.size());
    assert(i_maxx < points.size());
    assert(i_maxy < points.size());

    std::vector<Vector<2>> rtn;
    rtn.push_back(points[i_minx]);
    if (points[i_miny] != rtn[0])
        rtn.push_back(points[i_miny]);
    if ((points[i_maxx] != rtn[0]) && (points[i_maxx] != rtn[1]))
        rtn.push_back(points[i_maxx]);
    if ((points[i_maxy] != rtn[0]) && (points[i_maxy] != rtn[1]) && (points[i_maxy] != rtn[2]))
        rtn.push_back(points[i_maxy]);

    uint64_t num_cardinal_points = rtn.size();
    assert(0 < num_cardinal_points && num_cardinal_points <= 4);
    // std::cerr << "num cardinal points: " << num_cardinal_points << std::endl;
    
    switch (num_cardinal_points)
    {
    case 1:
        break;
    case 2:
        break;
    case 3:
        if (!MiscellaneousMath::is_convex_positively_oriented_polygon(rtn))
        {
            std::cerr << "not positively oriented triangle: " << std::endl;
            std::cerr << rtn[0].transpose() << std::endl;
            std::cerr << rtn[1].transpose() << std::endl;
            std::cerr << rtn[2].transpose() << std::endl;

            if (MiscellaneousMath::signed_area_of_triangle(std::array<Vector<2>, 3>({rtn[0], rtn[1], rtn[2]})) < 0)
                std::cerr << "wrong orientation" << std::endl;
            else if (MiscellaneousMath::signed_area_of_triangle(std::array<Vector<2>, 3>({rtn[0], rtn[1], rtn[2]})) == 0)
                std::cerr << "degenerate" << std::endl;
            else
                std::cerr << "WTF?!?!" << std::endl;
            assert(false);
        }
        assert(MiscellaneousMath::is_convex_positively_oriented_polygon(rtn));
        break;
    case 4:
        if (!MiscellaneousMath::is_convex_positively_oriented_polygon(rtn))
        {
            std::cerr << "not positively oriented quadrilateral: " << std::endl;
            std::cerr << rtn[0].transpose() << std::endl;
            std::cerr << rtn[1].transpose() << std::endl;
            std::cerr << rtn[2].transpose() << std::endl;
            std::cerr << rtn[3].transpose() << std::endl;
        }
        assert(MiscellaneousMath::is_convex_positively_oriented_polygon(rtn));
        break;
    
    default:
        assert(false); // should never get here
        break;
    }

    return rtn;
}

double MiscellaneousMath::signed_area_of_triangle(const std::array<Vector<2>,3>& triangle)
{
    Matrix<2> A;

    A.col(0) = triangle[1]-triangle[0];
    A.col(1) = triangle[2]-triangle[1];

    return A.determinant()/2;
}

bool MiscellaneousMath::is_convex_positively_oriented_polygon(const std::vector<Vector<2>>& polygon)
{
    if (polygon.size() <= 2)
        return false;

    std::array<Vector<2>,3> triangle;
    triangle[0] = polygon[0];
    triangle[1] = polygon[1];

    for (uint64_t i=2; i<polygon.size(); ++i)
    {
        triangle[i%3] = polygon[i];
        if (signed_area_of_triangle(triangle) <= 0)
            return false;
    }

    return true;
}

bool MiscellaneousMath::point_inside_convex_polygon(const Vector<2> p, const std::vector<Vector<2>>& convex_polygon)
{

    // std::cerr << "point: " << std::endl;
    // std::cerr << p.transpose() << std::endl;

    // std::cerr << "convex_polygon: " << std::endl;

    // for (auto p : convex_polygon)
    //     std::cerr << p.transpose() << std::endl;

    assert(MiscellaneousMath::is_convex_positively_oriented_polygon(convex_polygon));

    // std::cerr << "is convex_polygon ... OK " << std::endl;


    std::array<Vector<2>,3> triangle;
    triangle[0] = p;

    for (uint64_t i=0; i<convex_polygon.size(); ++i)
    {
        triangle[1] = convex_polygon[i];
        triangle[2] = convex_polygon[(i+1)%convex_polygon.size()];
        double signed_area = MiscellaneousMath::signed_area_of_triangle(triangle);

        // std::cerr << "triangle: "  << std::endl;
        // std::cerr << triangle[0].transpose() << triangle[1].transpose() << triangle[2].transpose() << std::endl;
        // std::cerr << "signed area: " << signed_area << std::endl;
        if (signed_area <= 0)
            return false;
    }

    return true;
}

std::vector<Vector<2>> MiscellaneousMath::convex_hull(const std::vector<Vector<2>>& points)
{
    assert(!points.empty());

    std::vector<Vector<2>> cardinal_points = distinct_cardinal_points(points);

    assert(!cardinal_points.empty());

    const Vector<2> anchor = cardinal_points.front(); // certainly on the convex hull
    // std::cerr << "anchor: " << anchor.transpose() << std::endl;

    std::vector<Vector<2>> path; 

    if (cardinal_points.size() >= 3)
    {
        // to improve performance, reduce number of points under consideration
        assert(is_convex_positively_oriented_polygon(cardinal_points));
        assert(!point_inside_convex_polygon(anchor, cardinal_points));

        for (auto p : points)
            if (!point_inside_convex_polygon(p, cardinal_points))
                path.push_back(p);
    }
    else
    {
        path = points;
    }

    if (path.size() < 2)
        return path;

    assert(path.size() >= 3);


    // std::cerr << "cardinal points: " << std::endl;

    // for (auto p : cardinal_points)
    //     std::cerr << p.transpose() << std::endl;


    // std::cerr << "path: " << std::endl;

    // for (auto p : path)
    //     std::cerr << p.transpose() << std::endl;


    // sort the path points in increasing order around anchor point

    struct CompareByAngleAroundAnchor
    {
        CompareByAngleAroundAnchor(const Vector<2>& anchor_point_) : anchor_point(anchor_point_){}

        bool operator()(const Vector<2>& v1, const Vector<2>& v2) const
        {
            const double tiny = std::numeric_limits<double>::min();
            if (v1.isApprox(v2, tiny))
                return false;
            if (v1.isApprox(anchor_point, tiny))
                return true;
            if (v2.isApprox(anchor_point, tiny))
                return false;

            // anchor_point, v1, v2 are distinct

            std::array<Vector<2>,3> triangle;
            triangle[0] = anchor_point;
            triangle[1] = v1;
            triangle[2] = v2;

            double signed_area = MiscellaneousMath::signed_area_of_triangle(triangle);

            if (signed_area > 0)
                return true;
            else if (signed_area < 0)
                return false;

            assert(signed_area == 0);
            // anchor_point, v1, v2 are collinear

            // since anchor point is assumed to be on the convex hull, anchor_point will not lie between v1 and v2
            return (anchor_point-v1).norm() < (anchor_point-v2).norm();
        }

        const Vector<2> anchor_point;
    };


    // CompareByAngleAroundAnchor comparer(anchor);
    // for (auto p : path)
    //     if (!(comparer(anchor, p) && !comparer(p, anchor)))
    //     {
    //         std::cerr << "nonminimality witness: " <<  p.transpose() << std::endl;
    //         assert(false);
    //     }

  std::sort(path.begin(), path.end(), CompareByAngleAroundAnchor(anchor));

//   std::sort(path.begin(), path.end(), comparer);
//    std::stable_sort(path.begin(), path.end(), comparer);

    // for (uint64_t i=1; i<path.size(); ++i)
    // {
    //     auto p = path[i-1];
    //     auto q = path[i];
    //     if (!comparer(p,q))
    //     {
    //         std::cerr << "missort at i = " << i << std::endl;
    //         std::cerr << "anchor: " << anchor.transpose() << std::endl;
    //         std::cerr << p.transpose() << std::endl;
    //         std::cerr << q.transpose() << std::endl;
    //         // if (p==q)
    //         //     std::cerr << "p==q" << std::endl;
    //         // if (p==anchor)
    //         //     std::cerr << "p==anchor" << std::endl;
    //         // if (anchor==q)
    //         //     std::cerr << "anchor==q" << std::endl;
                
    //         if (comparer(p, q))
    //             std::cerr << "p < q" << std::endl;
    //         if (comparer(q, p))
    //             std::cerr << "q < p" << std::endl;

    //         std::array<Vector<2>,3> triangle;
    //         triangle[0] = anchor;
    //         triangle[1] = p;
    //         triangle[2] = q;

    //         double signed_area = MiscellaneousMath::signed_area_of_triangle(triangle);
    //         std::cerr << "signed_area: " << signed_area << std::endl;


    //         assert(false);
    //     }
    // }



    // assert(comparer(anchor, path.front()));
    // assert(!comparer(path.front(), anchor));
    // assert(anchor.isApprox(path.front(), 1e-5));
    assert(anchor.isApprox(path.front(), std::numeric_limits<double>::min()));


    // if (!path.front().isApprox(anchor, std::numeric_limits<double>::min()))
    // {
    //     std::cerr << "anchor: " << anchor.transpose() << std::endl;
    //     std::cerr << "path.front(): " << path.front().transpose() << std::endl;


    //     for (auto p : path)
    //         std::cerr << p.transpose() << std::endl;



    //     assert(false);
    // }
    // else
    // {
    //     std::cerr << "OK!!!!";
    //     assert(false);
    // }

    std::stack<Vector<2>> hull;
    for (auto p : path)
    {
        hull.push(p);
        while (hull.size() >= 3 && !top_three_make_left_turn(hull))
            delete_one_below_top(hull);
    }

    // copy the stack contents to a vector (seems unnecessarily inefficient)
    std::vector<Vector<2>> rtn;
    while(!hull.empty())
    {
        rtn.push_back(hull.top());
        hull.pop();
    }
    std::reverse(rtn.begin(), rtn.end());

    assert(MiscellaneousMath::is_convex_positively_oriented_polygon(rtn));

    return rtn;
}

double MiscellaneousMath::area_of_polygon(const std::vector<Vector<2>>& polygon)
{
    std::array<Vector<2>,3> triangle;
    triangle[0] = polygon[0];

    double area = 0;
    for (uint64_t i=2; i<polygon.size(); ++i)
    {
        triangle[1] = polygon[i-1];
        triangle[2] = polygon[i];
        area += MiscellaneousMath::signed_area_of_triangle(triangle);
    }

    return area;
}

std::array<double, 3> MiscellaneousMath::min_max_variance(const std::vector<Vector<2>>& points)
{
    Vector<2> centroid = MiscellaneousMath::centroid<2>(points);
    std::vector<double> distances_to_centroid;
    for (auto point : points)
        distances_to_centroid.push_back((point-centroid).norm());

    double min_distance_to_centroid = *std::min_element(distances_to_centroid.begin(), distances_to_centroid.end());
    double max_distance_to_centroid = *std::max_element(distances_to_centroid.begin(), distances_to_centroid.end());

    assert(distances_to_centroid.size() >= 2);
    double variance = MiscellaneousMath::sample_variance(distances_to_centroid);

    return std::array<double, 3>({min_distance_to_centroid, max_distance_to_centroid, variance});
}

std::array<double, 2> MiscellaneousMath::tilted_bbox_dims(const std::vector<Vector<2>>& points)
{
    assert(!points.empty());
    const uint64_t n = points.size();

    Eigen::MatrixXd X(n,2);    
    for (uint64_t i=0; i<n; ++i)
        X.row(i) = points[i];

    auto mu = X.colwise().sum() / n;
    auto h = Eigen::MatrixXd::Ones(n, 1);    
    auto B = X - h*mu;

    assert(n > 1);
    auto C = B.transpose()*B/(n-1);

    Eigen::SelfAdjointEigenSolver<Matrix<2>> eigensolver(C);
    if (eigensolver.info() != Eigen::Success)
    {
        std::cerr << "eigensolver failed" << std::endl;
        assert(false);
    }

    Matrix<2> O;
    O.col(0) = eigensolver.eigenvectors().col(0);
    O.col(1) = eigensolver.eigenvectors().col(1);

    assert(O.isUnitary());

    std::vector<Vector<2>> transformed_points;
    for (auto point : points)
        transformed_points.push_back(O*point);

    auto bbox = MiscellaneousMath::bounding_box<2>(transformed_points);
    Vector<2> corner_min = bbox[0];
    Vector<2> corner_max = bbox[1];
    double x_width = corner_max(0) - corner_min(0);
    double y_width = corner_max(1) - corner_min(1);

    assert(x_width > 0);
    assert(y_width > 0);
    return std::array<double, 2>({x_width, y_width});
}
