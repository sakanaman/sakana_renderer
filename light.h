#ifndef LIGHT_H
#define LIGHT_H
#include "vec3.h"
#include <cmath>

class Light 
{
public:
    double i;
    Light(double i):i(i){};
    virtual Vec3 radi_intencity(double cos) = 0;
};

class Isotropic_Light : public Light 
{
public:
    Isotropic_Light(double i):Light(i){};
    Vec3 radi_intencity(double cos) 
    {
        return Vec3(i,i,i);
    }    
};

class Spot_Light : public Light
{
public:
    double s;
    Spot_Light(double i, double s):Light(i),s(s){};
    Vec3 radi_intencity(double cos) {
        return i * std::pow(cos,s) * Vec3(1.0, 1.0, 1.0);
    }
};

#endif