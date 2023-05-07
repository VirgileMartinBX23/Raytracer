#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <random> 

#define PI 3.141592653589793238


class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
//usefull functions
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double carre(double x){
    return x*x;
}

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

Vector random_cos(const Vector& N){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2*PI*r1)*sqrt(1-r2);
    double y = sin(2*PI*r1)*sqrt(1-r2);
    double z = sqrt(r2);

    Vector T1;
    if (N[0] == std::min(N[0], std::min(N[1], N[2]))){
        T1 = Vector(0, N[2], -N[1]);
        T1.normalize();
    }
    else if (N[1] == std::min(N[0], std::min(N[1], N[2]))){
        T1 = Vector(N[2], 0, -N[0]);
        T1.normalize();
    }
    else {
        T1 = Vector(N[1], -N[0], 0);
        T1.normalize();
    }
    Vector T2 = cross(N, T1);
    Vector V = x*T1 + y*T2 + z*N;
    return V;
}

#endif
