#ifndef RAY_H
#define RAY_H
#include "vec3.h"
#include <float.h>

class Ray
{
  public:
    Vec3 origin;
    Vec3 direction;
    double org[3];
    double dir[3];
    Ray(const Vec3 &_origin, const Vec3 &_direction) : origin(_origin), direction(_direction)
    {
      org[0] = _origin.x; org[1] = _origin.y; org[2] = _origin.z;
      dir[0] = _direction.x; dir[1] = _direction.y; dir[2] = _direction.z;
    };
};

bool intersectAABBvsRay(double aabb[2][3], const Ray& ray)
{
  double t_max = DBL_MAX;
  double t_min = -1*DBL_MAX;

  for (int i = 0; i < 3; i++)
  {
    double t1 = (aabb[0][i] - ray.org[i])/ray.dir[i];
    double t2 = (aabb[1][i] - ray.org[i])/ray.dir[i];
    double t_near = std::min(t1, t2);
    double t_far = std::max(t1, t2);
    t_max = std::min(t_max, t_far);
    t_min = std::max(t_min, t_near);

    if (t_min > t_max) return false;
  }
  return true;
}

#endif
