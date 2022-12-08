#include "Force/ForceSpring.h"
#include "Body/BodyCuboid.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceSpring::ForceSpring(
    Simulation* sim,
    const BodyCuboid* cuboid1, const BodyCuboid* cuboid2,
    dtype contact1, dtype contact2, dtype k) : Force(sim), _cuboid1(cuboid1), _cuboid2(cuboid2) {
    _contact1 = contact1;
    _contact2 = contact2;
    _k = k;
}

void ForceSpring::set_k(dtype k) {
    _k = k;
}

bool ForceSpring::on_cuboid(const BodyCuboid* cuboid, const Vector3& xw) {
    Vector3 xl = cuboid->_E_i0.topLeftCorner(3, 3) * xw + cuboid->_E_i0.topRightCorner(3, 1);
    for (int i = 0;i < 3;i++) {
        dtype dist = fabs(fabs(xl(i)) - cuboid->_length(i) / 2.);
        if (dist < 1e-5) {
            return true;
        }
    }
    return false;
}

void ForceSpring::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    std::cout << "Hello ComputeForce" << std::endl;
}

void ForceSpring::computeForceWithDerivative(
    std::cout << "Hello ComputeForceWithDerivative" << std::endl;
}

}