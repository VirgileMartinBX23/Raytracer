#ifndef RAY_H
#define RAY_H


#include "Vector.h"

// O is the starting point of the ray
// u is the unit vecteur of the ray 
class Ray {
public:
    Ray(const Vector& O, const Vector& u, double& n1): O(O), u(u), n1(n1){};
    Vector O, u;
    double n1;
};

#endif