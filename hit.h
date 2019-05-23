#ifndef HIT_H
#define HIT_H

#include "vec3.h"

class Figure;

class Hit
{ //衝突に関する情報の格納
  public:
    double t;       //衝突点までの距離
    Vec3 hitPos;    //衝突位置
    Vec3 hitNormal; //衝突位置の法線
    double u;
    double v;
    const Figure *hitShape;

    Hit()
    {
        t = 1000000;
        hitShape = nullptr;
    };
    Hit(double t, const Vec3 &hitPos, const Vec3 &hitNormal)
        : t(t), hitPos(hitPos), hitNormal(hitNormal){};
};

#endif