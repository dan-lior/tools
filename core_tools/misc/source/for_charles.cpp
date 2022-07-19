#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>

typedef std::array<double, 2> My2DPoint;

// returns x_min, y_min, x_max, y_max
std::array<double, 4> bounding_box(const std::vector<My2DPoint>& points)
{
    assert(!points.empty());

    std::vector<double> xs;
    std::vector<double> ys;
    for (auto p : points)
    {
        xs.push_back(p[0]);
        ys.push_back(p[1]);
    }

    std::array<double, 4> rtn;

    rtn[0] = *std::min_element(xs.begin(), xs.end());       // x_min
    rtn[1] = *std::min_element(ys.begin(), ys.end());       // y_min
    rtn[2] = *std::max_element(xs.begin(), xs.end());       // x_max
    rtn[3] = *std::max_element(ys.begin(), ys.end());       // y_max

    return rtn;
} 

double myOrientation(const My2DPoint& q, const My2DPoint& p1, const My2DPoint& p2)
{
    return (q[0] - p1[0])*(p2[1] - p1[1]) - (q[1] - p1[1])*(p2[0]-p1[0]);
}


// returns true if ray in the positive x-direction, starting at q, properly intersects the line segment p1, p2
bool ray_intersects_edge(const My2DPoint& q, const My2DPoint& p1, const My2DPoint& p2)
{
    // corner case: line segment is horizontal (proper intersection impossible)
    if (p1[1] == p2[1]) 
        return false;

    // check if the ray lies above the line segment or below the line segment 
    if (q[1] <= p1[1] && q[1] <= p2[1]) 
        return false;
    if (q[1] >= p1[1] && q[1] >= p2[1]) 
        return false;

    // check if p is on the same side of the line p1 p2 as the infinite end of the ray
    const My2DPoint infinity({std::numeric_limits<double>::max(), 0});
    if (myOrientation(q,p1,p2) * myOrientation(infinity,p1,p2) >= 0)
        return false;

    return true;
}


// simple implementation
// assumes that polygon does not intersect itself
std::vector<bool> simplest_point_in_polygon(const std::vector<My2DPoint>& query_points, const std::vector<My2DPoint>& polygon)
{
    std::vector<bool> rtn(query_points.size());

    auto bbox = bounding_box(query_points);
    double x_min = bbox[0];
    double y_min = bbox[1];
    double x_max = bbox[2];
    double y_max = bbox[3];

    for (auto q : query_points)
    {
        if (q[0] <= x_min || q[1] <= y_min || q[0] >= x_max || q[1] >= y_max) // q is not strictly within bounding box
            rtn.push_back(false);
        else
        {
            int num_intersections = 0;
            for (int i=0; i<polygon.size(); ++i)
                if (ray_intersects_edge(q, polygon[i], polygon[(i+1)%polygon.size()]))
                    num_intersections++;
                
            rtn.push_back(num_intersections % 2 == 1); // odd number of proper intersections means that point is inside polygon
        }
    }

    return rtn;
}