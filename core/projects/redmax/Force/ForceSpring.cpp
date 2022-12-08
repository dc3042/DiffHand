#include "Force/ForceSpring.h"
#include "Body/BodyCuboid.h"
#include "CollisionDetection/CollisionDetection.h"
#include "CollisionDetection/Contact.h"
#include "Body/Body.h"
#include "Simulation.h"
#include "Joint/Joint.h"

namespace redmax {

ForceSpring::ForceSpring(
    Simulation* sim,
    const BodyCuboid* cuboid1, const BodyCuboid* cuboid2,
    dtype contact1, dtype contact2, dtype k, dtype l) : Force(sim), _cuboid1(cuboid1), _cuboid2(cuboid2) {
    _contact1 = contact1;
    _contact2 = contact2;
    _k = k;
    _l = l;
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

    Matrix4 E1 = _cuboid1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _cuboid2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 xl1 = _cuboid1->get_contact_points()[_contact1];
    Vector3 xl2 = _cuboid1->get_contact_points()[_contact2];

    Vector3 xw1 = R1 * xl1 + p1;
    Vector3 xw2 = R2 * xl2 + p2;

    MatrixX G1 = math::gamma(xl1);
    MatrixX G2 = math::gamma(xl2);

    Vector3 dx = xw2 - xw1;
    dtype l = dx.norm();
    Vector3 f = _k * dx;

    fm.segment(_cuboid1->_index[0], 6) += G1.transpose() * (R1.transpose() * f);
    fm.segment(_cuboid2->_index[0], 6) -= G2.transpose() * (R2.transpose() * f);

    std::cout << "xw1 " << xw1 << std::endl;
    std::cout << "xw2 " << xw2 << std::endl;
    std::cout << "length " << (xw1 - xw2).norm() << std::endl;
    std::cout << "_l " << _l << std::endl;
    std::cout << "force " << f << std::endl;

}

void ForceSpring::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose) {
    std::cout << "Hello ComputeForceWithDerivative" << std::endl;

    Matrix4 E1 = _cuboid1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _cuboid2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 xl1 = _cuboid1->get_contact_points()[_contact1];
    Vector3 xl2 = _cuboid1->get_contact_points()[_contact2];

    Vector3 xw1 = R1 * xl1 + p1;
    Vector3 xw2 = R2 * xl2 + p2;

    MatrixX G1 = math::gamma(xl1);
    MatrixX G2 = math::gamma(xl2);

    Matrix3 I = Matrix3::Identity();

    Vector3 dx = xw2 - xw1;
    dtype l = dx.norm();
    Vector3 f = _k * dx;

    fm.segment(_cuboid1->_index[0], 6) += G1.transpose() * (R1.transpose() * f);
    fm.segment(_cuboid2->_index[0], 6) -= G2.transpose() * (R2.transpose() * f);

    std::cout << "xw1 " << xw1 << std::endl;
    std::cout << "xw2 " << xw2 << std::endl;
    std::cout << "length " << (xw1 - xw2).norm() << std::endl;
    std::cout << "_l " << _l << std::endl;
    std::cout << "force " << f << std::endl;

    dtype coeff = 1 - _l/l;

    Km.block(_cuboid1->_index[0], _cuboid1->_index[0], 6, 6) += 
            _k * G1.transpose() * (Matrix<dtype, 3, 6>() << coeff * math::skew(R1.transpose() * (xw2 - p1)), -I).finished();
    Km.block(_cuboid2->_index[0], _cuboid2->_index[0], 6, 6) +=
        _k * G2.transpose() * (Matrix<dtype, 3, 6>() << coeff * math::skew(R2.transpose() * (xw1 - p2)), -I).finished();
    Km.block(_cuboid1->_index[0], _cuboid2->_index[0], 6, 6) +=
        _k * G1.transpose() * (R1.transpose() * (R2 * (Matrix<dtype, 3, 6>() << -coeff * math::skew(xl2), I).finished()));
    Km.block(_cuboid2->_index[0], _cuboid1->_index[0], 6, 6) += 
        _k * G2.transpose() * (R2.transpose() * (R1 * (Matrix<dtype, 3, 6>() << -coeff * math::skew(xl1), I).finished()));
}

}