#ifndef FIGURE_H
#define FIGURE_H

#include "vec3.h"
#include "texture.h"
#include "ray.h"
#include "hit.h"
#include "light.h"
#include <cmath>
#include <memory>
#include <utility>

class Figure
{
  public:
    int material;
    std::shared_ptr<Texture> texture;
    double bbox[2][3];
    std::shared_ptr<Light> light; 
    Figure(){};
    Figure(std::shared_ptr<Texture> texture, int material, std::shared_ptr<Light> light) : texture(texture), material(material), light(light){};
    virtual bool intersect(const Ray &ray, Hit &hit) const = 0;
    virtual double area() const = 0;
    virtual Vec3 randompoint() const = 0; 
};

class Sphere : public Figure
{
  public:
    Vec3 center;
    double radius;
    Sphere(const Vec3 &center, double radius, std::shared_ptr<Texture> texture, int material, std::shared_ptr<Light> light)
        : center(center), radius(radius), Figure(texture, material, light)
        {
            bbox[0][0] = center.x - radius;
            bbox[1][0] = center.x + radius;
            bbox[0][1] = center.y - radius;
            bbox[1][1] = center.y + radius;
            bbox[0][2] = center.z - radius;
            bbox[1][2] = center.z + radius;
        };

    virtual bool intersect(const Ray &ray, Hit &hit) const
    {
        double d_norm = ray.direction.length();
        double oc_norm = (ray.origin - center).length();

        double a = d_norm * d_norm;
        double b = 2 * dot(ray.direction, ray.origin - center);
        double c = oc_norm * oc_norm - radius * radius;
        double d = b * b - 4 * a * c;
        if (d < 0)
            return false;

        double t1 = (-b - sqrt(d)) / (2 * a);
        double t2 = (-b + sqrt(d)) / (2 * a);

        double t = t1;
        if (t <= 1e-6)
        {
            t = t2;
            if (t <= 1e-6)
                return false;
        }

        hit.t = t;
        hit.hitPos = ray.origin + t * ray.direction;
        hit.hitNormal = normalize(hit.hitPos - center);
        hit.hitShape = this;
        double phi = std::atan2(hit.hitNormal.z, hit.hitNormal.x) - M_PI / 2.0;
        if (phi < 0)
            phi += 2 * M_PI;
        if (phi > 2 * M_PI)
            phi -= 2 * M_PI;
        double theta = std::acos(hit.hitNormal.y);
        hit.u = phi / (2 * M_PI);
        hit.v = theta / M_PI;
        return true;
    };

    virtual double area() const 
    {
        return 4*M_PI*radius*radius;
    }

    virtual Vec3 randompoint() const
    {
        return center + radius * randomsphere();
    }
};

void renritsu_3(const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &origin, const Vec3 &direction, double &m, double &n, double &t)
{
    Vec3 E1 = b - a;
    Vec3 E2 = c - a;
    Vec3 T = origin - a;
    Vec3 R = direction;

    Vec3 P = cross(R, E2);
    Vec3 Q = cross(T, E1);

    t = dot(Q, E2) / dot(P, E1);
    m = dot(P, T) / dot(P, E1);
    n = dot(Q, R) / dot(P, E1);
}

class Tri : public Figure
{
  public:
    Vec3 a;
    Vec3 b;
    Vec3 c;

    Tri(const Vec3 &a, const Vec3 &b, const Vec3 &c, std::shared_ptr<Texture> texture, int material, std::shared_ptr<Light> light) : a(a), b(b), c(c), Figure(texture, material, light)
    {
        bbox[1][0] = std::max(a.x,std::max(b.x, c.x));
        bbox[0][0] = std::min(a.x,std::min(b.x, c.x));
        bbox[1][1] = std::max(a.y, std::max(b.y, c.y));
        bbox[0][1] = std::min(a.y, std::min(b.y,c.y));
        bbox[1][2] = std::max(a.z,std::max(b.z, c.z));
        bbox[0][2] = std::min(a.z, std::min(b.z, c.z));
    };

    virtual bool intersect(const Ray &ray, Hit &hit) const
    {
        double m, n, t;
        renritsu_3(a, b, c, ray.origin, ray.direction, m, n, t);
        if (m + n <= 1 && 0 <= m && m <= 1 && n <= 1 && n >= 0 && t > 1e-6)
        {
            hit.t = t;
            hit.hitPos = ray.origin + t * ray.direction;
            hit.hitNormal = -1*normalize(cross(b - a, c - a));

            hit.hitShape = this;
            hit.u = n / (m + n);
            hit.v = m + n;

            return true;
        }
        return false;
    }

    virtual double area() const
    {
        return 0.5 * ( cross(b-a,c-a).length() );
    }

    virtual Vec3 randompoint() const
    {
        return randomtriangle(a,b,c);
    }
};


#endif
