#pragma once

#include "common_defs.h"
#include "misc.h"
#include "misc_math.h"
#include "affinity.h"

struct Polyline
{
    Polyline(const std::vector<Vector<3>>& vertices, uint64_t smoothing_radius_=1);
    
    Polyline(const std::string& filename_vtk, uint64_t smoothing_radius_=1);

    Polyline reverse_direction() const;

    std::vector<double> distances_from_start() const;

    std::vector<double> curvature() const;

    // current implementation recomputes curvatures
    std::vector<double> torsion() const;

    // resample uniformly so that mm is the approximate distance between samples
    Polyline resample(double mm) const;

    Polyline reset_smoothing_radius(uint64_t new_smoothing_radius) const;

    Polyline transform(const AffinityIso<3>& affinity) const;

    Polyline gaussian_perturbation(double mu, double sigma) const;

    void display_stats() const;

    void write_to_file(const std::string& filename) const;

    static AffinityIso<3> register_icp(const Polyline& fixed, const Polyline& floating, double convergence_tolerance, uint64_t max_iterations, bool vanilla);

    std::vector<Vector<3>> vertices;
    std::vector<Vector<3>> tangent_directions;

    private:

//    const uint64_t smoothing_radius; // compiler error?
    uint64_t smoothing_radius;

    void generate_smooth_tangents();

};
