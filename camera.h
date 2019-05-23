#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"
#include "vec3.h"
#include "myrandom.h"

class Camera
{
  public:
    Vec3 camPos;
    Vec3 camForward;
    Vec3 camRight;
    Vec3 camUp;

    Camera(const Vec3 &_camPos, const Vec3 &_camForward) : camPos(_camPos), camForward(_camForward)
    {
        camRight = -1 * normalize(cross(camForward, Vec3(0, 1, 0)));
        camUp = normalize(cross(camForward, camRight));
    };

    virtual Ray getRay(double u, double v, double& w) const = 0;
};

class Pinhole_Camera : public Camera
{
  public:
    Pinhole_Camera(const Vec3 &_camPos, const Vec3 &_camForward) : Camera(_camPos, _camForward){};

    Ray getRay(double u, double v, double& w) const
    {
        u = -u;
        Vec3 pinhole = camPos + camForward;
        Vec3 sensorPos = camPos + u * camRight + v * camUp;
        w = pow( std::abs(dot(normalize(pinhole - sensorPos),normalize(camForward))) ,4);
        return Ray(sensorPos, normalize(pinhole - sensorPos));
    };
};

class Thin_Lens_Camera : public Camera
{
  public:
    double a;
    double b;
    double f;
    double R;
    double F_Number;
    Vec3 center;
    Vec3 forcus;
    Thin_Lens_Camera(const Vec3 &camPos, const Vec3 &camForward, double a, const Vec3 &forcus, double F_Number) : Camera(camPos, camForward), a(a), forcus(forcus), F_Number(F_Number)
    {
        center = camPos + normalize(camForward) * a;
        b = (forcus - camPos).length() - a;
        f = (a * b) / (a + b);
        R = f / (2.0 * F_Number);
    };

    Ray getRay(double u, double v, double& w) const
    {
        u = -u;
        Vec3 sensorPos = camPos + u * camRight + v * camUp;
        Vec3 L = center + R * random_unit_circle(normalize(camForward));
        Vec3 l = normalize(center - sensorPos);
        Vec3 P = sensorPos + ((a + b) / dot(normalize(camForward), l)) * l;
        w = 1;
        return Ray(L, normalize(P - L));
    }
};

#endif