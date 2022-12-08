#pragma once
#include "Common.h"
#include "Utils.h"
#include "Force/Force.h"

namespace redmax {

class Body;

class ForceSpring : public Force {
public:
    const BodyCuboid* _cuboid1; // cuboid bodies to consider contact
    const BodyCuboid* _cuboid2;
    dtype _contact1; // cuboid contact points
    dtype _contact2;
    dtype _k;              // spring constant
    

    ForceSpring(
        Simulation* sim,
        const BodyCuboid* cuboid1, const BodyCuboid* cuboid2, dtype contact1, dtype contact2, dtype k = 1.);

    void set_k(dtype _k);

    bool on_cuboid(const BodyCuboid* cuboid, const Vector3& xw);

    void computeForce(VectorX& fm, VectorX& fr, bool verbose = false);
    void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose = false);
};

}