#include "quadric.h"


std::vector<Quadric> Quadric::fit_all_types_of_quadrics(const std::string& filename)
{
    std::cerr << "fitting from " << filename << std::endl;

    // all_quadrics::TriangleMesh inputMesh;

    // if (!inputMesh.loadObj(filename.c_str())) {
    //     std::cerr << "Couldn't load file " << filename << endl;
    //     assert(false);
    // }

    // //Always recenter and scale your data before fitting!
    // // inputMesh.centerAndScale(1);
    // auto adjustment = inputMesh.dan_centerAndScale(1);
    // vec3 center = adjustment.first;
    // double scale = adjustment.second;

    // std::vector<quadric_fit::Quadric> qfits;
    // quadric_fit::fitAllQuadricTypes(inputMesh, qfits);
    // assert(qfits.size() == QuadricFitTypes::NUM_QUADRIC_TYPES);
    

    std::vector<Quadric> rtn;
    // for (auto qfit : qfits)
    // {
    //     qfit.normalizeField(); // is this necessary?
    //     auto q = qfit.q;
    //     rtn.push_back(Quadric(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9]));
    // }

    return rtn;
}
