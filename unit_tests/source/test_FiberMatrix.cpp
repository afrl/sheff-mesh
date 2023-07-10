#include "mesh/algorithms/FiberMatrix.h"

using namespace mesh;

int main()
{
    // Create RVE that looks like:
    //    _
    // o-/-\----------o
    // ||   |         |
    // | \_/     _    |
    // |        / \   |
    // |       |   |  |
    // |        \_/   |
    // |             /|\
    // o------------|-o |
    //               \_/



    FiberMatrixRVEDescription rve;
    rve.fibers = {
        { { 0.3, 0.8, 0 }, 0.25 },
        { { 0.7, 0.45, 0 }, 0.25 },
        { { 1, 0, 0 }, 0.25 },
    };
    rve.rve_boundary = make_rve_boundary({ { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } }, 0.0);
    mesh_rve(rve, 0.1);

    return 0;
}