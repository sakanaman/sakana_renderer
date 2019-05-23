#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <cmath>
using namespace std;

class Vec3
{
  public:
    double x, y, z;
    Vec3()
    {
        x = 0;
        y = 0;
        z = 0;
    };
    Vec3(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    };
    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    };
    double length2() const
    {
        return x * x + y * y + z * z;
    };
};

Vec3 operator+(const Vec3 &left, const Vec3 &right)
{
    return Vec3(left.x + right.x, left.y + right.y, left.z + right.z);
}
Vec3 operator-(const Vec3 &left, const Vec3 &right)
{
    return Vec3(left.x - right.x, left.y - right.y, left.z - right.z);
}
Vec3 operator*(double k, const Vec3 &v)
{
    return Vec3(k * v.x, k * v.y, k * v.z);
}
Vec3 operator*(const Vec3 &v, double k)
{
    return k * v;
}
Vec3 operator*(const Vec3 &v1, const Vec3 &v2)
{
    return Vec3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

Vec3 operator/(double k, const Vec3 &v)
{
    return Vec3(k / v.x, k / v.y, k / v.z);
}
Vec3 operator/(const Vec3 &v, double k)
{
    return Vec3(v.x / k, v.y / k, v.z / k);
}
bool operator==(const Vec3 &v1, const Vec3& v2)
{
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

bool operator!=(const Vec3& v1, const Vec3& v2)
{
    return !(v1==v2);
}

std::ostream &operator<<(std::ostream &stream, const Vec3 &v)
{
    stream << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return stream;
}
double dot(const Vec3 &v1, const Vec3 &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
Vec3 cross(const Vec3 &v1, const Vec3 &v2)
{
    return Vec3(v1.y * v2.z - v2.y * v1.z, v1.z * v2.x - v2.z * v1.x, v1.x * v2.y - v2.x * v1.y);
}

Vec3 normalize(const Vec3 &v)
{
    return v / v.length();
}
Vec3 reflect(const Vec3 &d, const Vec3 &n)
{
    return d - 2 * dot(d, n) * n;
}

void orthonormalBasis(const Vec3 &n, Vec3 &x, Vec3 &z)
{
    if (n.x > 0.9)
        x = Vec3(0, 1, 0);
    else
        x = Vec3(1, 0, 0);
    x = normalize(x - dot(x, n) * n);
    z = normalize(cross(n, x));
}
void coordinate_system(const Vec3 &v1, Vec3 &v2, Vec3 &v3)
{
    if (std::abs(v1.x) > std::abs(v1.y))
        v2 = Vec3(-v1.z, 0, v1.x) /
                std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        v2 = Vec3(0, v1.z, -v1.y) /
             std::sqrt(v1.y * v1.y + v1.z * v1.z);
    v3 = normalize(cross(v2, v1));//v2とv1入れ替えた.
}
void Transpose3x3(Vec3 &a, Vec3 &b, Vec3 &c)
{
    double p = a.y;
    double q = a.z;
    double r = b.z;
    a.y = b.x;
    a.z = c.x;
    b.z = c.y;
    b.x = p;
    c.x = q;
    c.y = r;
}

inline bool isNan(const Vec3 &v)
{
    return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z);
}

#endif
