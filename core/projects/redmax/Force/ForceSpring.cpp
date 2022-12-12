#include "Force/ForceSpring.h"
#include "Body/Body.h"
#include "Simulation.h"

namespace redmax {

ForceSpring::ForceSpring(
    Simulation* sim,
    const Body* body1, const Body* body2,
    dtype contact1, dtype contact2, dtype k, dtype l) : Force(sim), body1(body1), body2(body2) {
    _contact1 = contact1;
    _contact2 = contact2;
    _k = k;
    _l = l;
}

void ForceSpring::set_k(dtype k) {
    _k = k;
}


void ForceSpring::computeForce(VectorX& fm, VectorX& fr, bool verbose) {
    //std::cout << "Hello ComputeForce" << std::endl;

    Matrix4 E1 = _body1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _body2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 xl1 = _body1->get_contact_points()[_contact1];
    Vector3 xl2 = _body1->get_contact_points()[_contact2];

    Vector3 xw1 = R1 * xl1 + p1;
    Vector3 xw2 = R2 * xl2 + p2;

    MatrixX G1 = math::gamma(xl1);
    MatrixX G2 = math::gamma(xl2);

    Vector3 dx = xw2 - xw1;
    dtype l = dx.norm();

    dtype f =  - _k * (l - _l);


    fm.segment(_body1->_index[0], 6) +=  f/l * G1.transpose() * (R1.transpose() * dx);
    fm.segment(_body2->_index[0], 6) -=  f/l * G2.transpose() * (R2.transpose() * dx);

    //std::cout << "xw1 " << xw1 << std::endl;
    //std::cout << "xw2 " << xw2 << std::endl;
    //std::cout << "length " << l << std::endl;
    //std::cout << "_l " << _l << std::endl;
    //std::cout << "force direction" << f.dot(dx) << std::endl;
    //exit(0);

}

void ForceSpring::computeForceWithDerivative(
    VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose) {
    //std::cout << "Hello ComputeForceWithDerivative" << std::endl;

    Matrix4 E1 = _body1->_E_0i;
    Matrix3 R1 = E1.topLeftCorner(3, 3);
    Vector3 p1 = E1.topRightCorner(3, 1);
    Matrix4 E2 = _body2->_E_0i;
    Matrix3 R2 = E2.topLeftCorner(3, 3);
    Vector3 p2 = E2.topRightCorner(3, 1);

    Vector3 xl1 = _body1->get_contact_points()[_contact1];
    Vector3 xl2 = _body1->get_contact_points()[_contact2];

    Vector3 xw1 = R1 * xl1 + p1;
    Vector3 xw2 = R2 * xl2 + p2;

    MatrixX G1 = math::gamma(xl1);
    MatrixX G2 = math::gamma(xl2);

    Matrix3 I = Matrix3::Identity();

    Vector3 dx = xw2 - xw1;
    dtype l = dx.norm();

    dtype f =  - _k * (l - _l); // No dependence on l_dot

    fm.segment(_body1->_index[0], 6) +=  f/l * G1.transpose() * R1.transpose() * dx;
    fm.segment(_body2->_index[0], 6) -=  f/l * G2.transpose() * R2.transpose() * dx;

    std::cout << "xw1 " << xw1 << std::endl;
    std::cout << "xw2 " << xw2 << std::endl;
    std::cout << "length " << l << std::endl;
    std::cout << "_l " << _l << std::endl;
    std::cout << "force " << f << std::endl;

    /**
    K1
    */

    dtype df_dl = -_k;

    VectorX dl_dq = VectorX::Zero(12);
    dl_dq.segment(0, 6) = - dx.transpose()/ l * R1 * G1;
    dl_dq.segment(6, 6) = dx.transpose() / l * R2 * G2;
    VectorX df_dq = df_dl * dl_dq;
    VectorX dfl_dq =  1/l * df_dq - f/(l*l) * dl_dq;

    VectorX f_vector = VectorX::Zero(12);
    f_vector.segment(0, 6) = G1.transpose() * R1.transpose() * dx;
    f_vector.segment(6, 6) = - G2.transpose() * R1.transpose() * dx;

    MatrixX K1 = dfl_dq * f_vector.transpose();

    Km.block(_body1->_index[0], _body1->_index[0], 6, 6) += K1.block(0, 0, 6, 6);
    Km.block(_body1->_index[0], _body2->_index[0], 6, 6) += K1.block(0, 6, 6, 6);
    Km.block(_body2->_index[0], _body1->_index[0], 6, 6) += K1.block(6, 0, 6, 6);
    Km.block(_body2->_index[0], _body2->_index[0], 6, 6) += K1.block(6, 6, 6, 6);

    /**
    K2
    */

    Km.block(_body1->_index[0], _body1->_index[0], 3,3) += f/l * math::skew(xl1) * math::skew(R1.transpose() * (p1 - xw2));
    Km.block(_body1->_index[0], _body1->_index[3], 3,3) += f/l * math::skew(xl1);
    Km.block(_body1->_index[0], _body2->_index[0], 3,3) += f/l * math::skew(xl1) * R1.transpose() * R2 * math::skew(xl2);
    Km.block(_body1->_index[0], _body2->_index[3], 3,3) += f/l * -math::skew(xl1) * R1.transpose() * R2;

    Km.block(_body1->_index[3], _body1->_index[0], 3,3) += f/l * math::skew(R1.transpose() * (p1 - xw2));
    Km.block(_body1->_index[3], _body1->_index[3], 3,3) += f/l * I;
    Km.block(_body1->_index[3], _body2->_index[0], 3,3) += f/l * R1.transpose() * R2 * math::skew(xl2);
    Km.block(_body1->_index[3], _body2->_index[3], 3,3) += f/l * -R1.transpose() * R2;

    Km.block(_body2->_index[0], _body1->_index[0], 3,3) += f/l * math::skew(xl2) * R2.transpose() * R1 * math::skew(xl1);
    Km.block(_body2->_index[0], _body1->_index[3], 3,3) += f/l * -math::skew(xl2) * R2.transpose() * R1;
    Km.block(_body2->_index[0], _body2->_index[0], 3,3) += f/l * math::skew(xl2) * math::skew(R2.transpose() * (p2 - xw1));
    Km.block(_body2->_index[0], _body2->_index[3], 3,3) += f/l * math::skew(xl2);

    Km.block(_body2->_index[3], _body1->_index[0], 3,3) += f/l * R2.transpose() * R1 * math::skew(xl1);
    Km.block(_body2->_index[3], _body1->_index[3], 3,3) += f/l * -R2.transpose() * R1;
    Km.block(_body2->_index[3], _body2->_index[0], 3,3) += f/l * math::skew(R2.transpose() * (p2 - xw1));
    Km.block(_body2->_index[3], _body2->_index[3], 3,3) += f/l * I;

}

}