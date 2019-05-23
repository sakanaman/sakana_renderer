#ifndef MYRANDOM_H
#define MYRANDOM_H

#include <random>
#include <cmath>
#include "vec3.h"

std::random_device dev_rnd;
std::mt19937 mt(dev_rnd());
std::uniform_real_distribution<> dist1(0, 1);
std::uniform_real_distribution<> dist2(-1,1);

inline double rnd()
{
    return dist1(mt);
}

inline double rnd2()
{
    return dist2(mt);
}

Vec3 random_unit_circle(const Vec3 &n)
{
    double u = rnd();
    double r = rnd();

    double y = 0;
    double x = r * std::cos(2 * M_PI * u);
    double z = r * std::sin(2 * M_PI * u);

    Vec3 xv, zv;
    coordinate_system(n, xv, zv);

    return x * xv + y * n + z * zv;
}

Vec3 randomHemisphere(const Vec3 &n)
{
    double u = rnd();
    double v = rnd();

    double y = u;
    double x = std::sqrt(1 - u * u) * std::cos(2 * M_PI * v);
    double z = std::sqrt(1 - u * u) * std::sin(2 * M_PI * v);

    Vec3 xv, zv;
    coordinate_system(n, xv, zv);

    return x * xv + y * n + z * zv; //{xv,n,v}を基底とした時の[x,y,z]
}

Vec3 randomsphere()
{
    double y = rnd2();
    double phi = rnd() * 2 * M_PI;
    double x = sqrt(1-y*y)*cos(phi);
    y = y;
    double z = sqrt(1-y*y)*sin(phi);
    return Vec3(x,y,z);
}

Vec3 randomtriangle(const Vec3& a, const Vec3& b, const Vec3& c)
{
    double u1 = rnd();
    double u2 = rnd();

    return a*(1-sqrt(u1)) + b*( sqrt(u1)*(1-u2) ) + c*sqrt(u1)*u2;      
}

Vec3 randomCosineHemisphere(const Vec3 &n)
{
    double r1 = rnd();
    double r2 = rnd();
    double phi = 2 * M_PI * r1;
    double y = sqrt(1-r2);
    double x = std::cos(phi)* sqrt(r2);
    double z = std::sin(phi)* sqrt(r2);
    Vec3 xv, zv;
    coordinate_system(n, xv, zv);
    return normalize(x * xv + y * n + z * zv);
}

#endif