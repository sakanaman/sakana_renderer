#ifndef GGX_H
#define GGX_H
#include "vec3.h"
double chi_plus(double x)
{
    if(x>0) return 1.0;
    return 0;
}

double F(const double ni, const Vec3& o, const Vec3& h)
{
    double c = dot(o,h);
    double g = ni*ni + c*c - 1;
    return 0.5 * (g-c)*(g-c) / (g+c)/(g+c) * 
            (1 + (c*(g+c)-1)*(c*(g+c)-1)/(c*(g-c)+1)/(c*(g-c)+1));
}

double F_die(const double& cos, const double n_ob, const double n_v) {
    double R = (n_ob - n_v)*(n_ob - n_v) / (n_ob + n_v)/(n_ob + n_v);
    return R + (1 - R) * std::pow(1-cos, 5);
}
Vec3 F(const double cost, const Vec3& R)
{
    return R + (Vec3(1,1,1) - R) * std::pow(1 - cost, 5);
}

double T( const double ni, const Vec3& i, const Vec3& h)
{
    return 1 - F(ni, i, h);
}

double D_GGX(const Vec3& m, const Vec3& n, double alpha) 
{
    return alpha*alpha/M_PI/std::abs(dot(m,n))/std::abs(dot(m,n))/std::abs(dot(m,n))/std::abs(dot(m,n))
            /(alpha*alpha + tan(acos(dot(m,n)))*tan(acos(dot(m,n)))) / (alpha*alpha + tan(acos(dot(m,n)))*tan(acos(dot(m,n))));
}

double G1_GGX(const Vec3& u, const Vec3& n, double alpha) 
{
    return 2/(1 + sqrt(1 + alpha*alpha*tan(acos(dot(u,n)))*tan(acos(dot(u,n))))); 
}

double G(const Vec3& i, const Vec3& o, const Vec3& n, double alpha)
{
    return G1_GGX(i,n,alpha) * G1_GGX(o,n,alpha);
}
#endif