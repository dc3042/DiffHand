#pragma once
#include "Common.h"
#include "Utils.h"
#include "Force/Force.h"
#include "Body/Body.h"

namespace redmax {

class Body;

class ForceSpring : public Force {
public:
    const Body* _body1; // cuboid bodies to consider contact
    const Body* _body2;
    dtype _contact1; // cuboid contact points
    dtype _contact2;
    dtype _k;              // spring constant
    dtype _l;              // spring default length
    

    ForceSpring(
        Simulation* sim,
        const Body* body1, const Body* body2, 
        dtype contact1, dtype contact2, 
        dtype k = 1., dtype l = 1.);

    void set_k(dtype _k);

    void computeForce(VectorX& fm, VectorX& fr, bool verbose = false);
    void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose = false);
};

}