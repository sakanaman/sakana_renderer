#ifndef RADI_H
#define RADI_H
#include <vector>
#include <algorithm>
#include "vec3.h"
#include "ibl.h"
#include "accel.h"
#include "ray.h"
#include "bvh.h"
#include "myrandom.h"
#include "hit.h"
#include "brdf.h"

#ifdef USE_DEBUG
Vec3 Debug_ao(const Ray& ray, Accel& accel) {
    Ray nextRay = ray;
    Vec3 I(0,0,0), weight(1,1,1);
    return I;
    // !implementing now!
}
Vec3 Debug_normal(const Ray& ray, Accel& accel) {
    Ray nextRay = ray;
    Vec3 I(0,0,0);

    std::shared_ptr<Figure> shape;
    Hit hit;
    shape = accel.intersect(nextRay, hit, 0, nodes);

    I.x = 0.5 + 0.5 * hit.hitNormal.x;
    I.y = 0.5 + 0.5 * hit.hitNormal.y;
    I.z = 0.5 + 0.5 * hit.hitNormal.z;

    return I;
}
Vec3 Debug_depth(const Ray& ray, Accel& accel) {
    Ray nextRay = ray;
    int depth = 0;
    double Prr = 1.0;
    Vec3 dir;
    std::shared_ptr<Figure> shape;
    for(;;depth++) {
        Hit hit;
        shape = accel.intersect(nextRay, hit, 0, nodes);
        if(shape != nullptr) {
            Vec3 orienting_normal = dot(nextRay.direction, hit.hitNormal) > 0 ? -1*hit.hitNormal : hit.hitNormal;
            dir = randomHemisphere(orienting_normal);
            nextRay = Ray(hit.hitPos, dir);
            double u = rnd();
            if(Prr < rnd()) {
                break;
            }
            Prr *= 0.8;
        }
        else{
            break;
        }
    }
    return Vec3(1,1,1) * Prr;
}
#endif

Vec3 getColor(const Ray &ray,const IBL &ibl,Accel& accels)
{
    Vec3 I = Vec3(0,0,0);
    double Prr = 1.0;
    Vec3 weight = Vec3(1.0,1.0,1.0);
    Ray nextRay = ray;
    std::shared_ptr<Figure> shape;
    double costerm=1.0,pdf=1.0;
    Vec3 bsdf,dir;
    bool is_specular = false;
    for(int depth=0;;depth++)
    {   
        Hit hit;
        shape =  accels.intersect(nextRay,hit,0,nodes);
        if( shape != nullptr )
        {
            if(hit.hitShape->light != nullptr && dot(hit.hitNormal,nextRay.direction) < 0)//光源にヒット
            {
                double c_o_s = std::abs(dot(hit.hitNormal, nextRay.direction));
                ////////////NEE
                if(is_specular == true)//完全スペキュラーをトレースしてきた
                {
                    I = I + weight * hit.hitShape->light->radi_intencity(c_o_s);
                    is_specular = false;
                    break;
                }
                if(depth == 0)//しょっぱなヒットしたx1が光源だった.
                {
                    I = I + weight * hit.hitShape->light->radi_intencity(c_o_s);
                    break;
                }
                /////////////
                double l = (hit.hitPos - nextRay.origin).length();
                double light_pdf = 1/accels.light_area * l * l / costerm;
                double mis_weight = pdf/(light_pdf + pdf);
                I = I + mis_weight * weight * hit.hitShape->light->radi_intencity(c_o_s);//NEEの時はコメントアウト
            }
            
            const Vec3 orienting_normal = dot(-1*nextRay.direction,hit.hitNormal)>0 ? hit.hitNormal : -1.0*hit.hitNormal;
            const Vec3 nowobjcolor = hit.hitShape->texture->getColor1(hit.u,hit.v,hit.hitPos);
            ////////マテリアルごとにbsdf,dir(入射光の方向),pdf,cos項を求める.
            switch(hit.hitShape->material)
            {
                case 0://diffuse
                {
                    bsdf = nowobjcolor * (1/M_PI);

                    do
                    {
                        dir = randomCosineHemisphere(orienting_normal);
                    }while(dot(dir,orienting_normal) <= 0);
                    costerm = std::abs(dot(dir,orienting_normal));
                    pdf = costerm/M_PI;
                    is_specular = false;

                    /////////NEEに関する処理
                    //double pa = 1/accels.light[i]->area();
                    int i = accels.select_light();
                    double pa =1/accels.light_area;
                    Vec3 x_l = accels.light[i]->randompoint();
                    Ray shadowRay = Ray(hit.hitPos, normalize(x_l - hit.hitPos));
                    double length_xl_x = (x_l-hit.hitPos).length();
                    Hit shadow_hit;
                    if(accels.intersect(shadowRay, shadow_hit, 0, nodes) != nullptr)
                    {
                        if( std::abs(shadow_hit.t - length_xl_x) < 1e-8 && dot(shadowRay.direction,shadow_hit.hitNormal) < 0 && dot(orienting_normal, shadowRay.direction) > 1e-8)
                        {   
                            double bsdf_pdf = dot(shadowRay.direction, orienting_normal)/M_PI * dot(shadowRay.direction, orienting_normal)/length_xl_x/length_xl_x;
                            double mis_weight = pa/(bsdf_pdf + pa);

                            double G = std::abs(dot(shadow_hit.hitNormal,normalize(hit.hitPos - x_l)))*std::abs(dot(orienting_normal,normalize(x_l - hit.hitPos))) / pow(length_xl_x,2);
                            I = I +  mis_weight * weight * accels.light[i]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * bsdf * G *(1/pa);
                        }
                    }
                    /////////
                    break;
                }
                case 1://specular
                    dir = reflect(nextRay.direction,orienting_normal);
                    costerm = std::abs(dot(dir,orienting_normal));
                    bsdf = nowobjcolor*(1/costerm);
                    pdf = 1.0;
                    is_specular = true;
                    break;

                case 2://dielectric
                {   
                    is_specular = true;
                    double u = rnd();
                    bool in = (dot(hit.hitNormal, orienting_normal) > 0);
                    double nt = in? 1.5 : 1.0;
                    double ni = in? 1.0 : 1.5;
                    double costhetaI = dot(orienting_normal, -1 * nextRay.direction);
                    double sinthetaT = ni/nt * std::sqrt(std::max(1 - costhetaI * costhetaI,0.0));
                    if(sinthetaT >= 1) { // total reflect
                        pdf = 1.0;
                        dir = reflect(nextRay.direction, orienting_normal);
                        //if(std::abs( dot(dir, orienting_normal) - dot(orienting_normal, -1*nextRay.direction)) >= 1e-8) std::cout << "fail reflect" << std::endl;
                        costerm = std::abs(dot(dir, orienting_normal));
                        bsdf = nowobjcolor / costerm;
                    }
                    else {
                        double costhetaT = std::sqrt(1 - sinthetaT*sinthetaT);
                        double Fr = F_die(costhetaI, 1.0 , 1.5);
                        double u = rnd();
                        if(u < Fr) {
                            pdf = Fr;
                            dir = reflect(nextRay.direction, orienting_normal);
                            costerm = std::abs(dot(dir, orienting_normal));
                            bsdf = Fr * nowobjcolor / costerm;
                        }
                        else {
                            pdf = 1.0 - Fr;
                            dir = refract(-1*nextRay.direction, orienting_normal, costhetaI, costhetaT, ni/nt);
                            costerm = std::abs(dot(dir, orienting_normal));
                            bsdf = (ni*ni)/(nt*nt) * (1.0 - Fr) / costerm * nowobjcolor;
                        }
                    }
                    break;
                }
                case 3://GGX - rough - specular
                {
                    double alpha_g = 0.4 * 0.4;
                    Vec3 u,v,halfv;
                    coordinate_system(orienting_normal, u, v);
                    double phi, theta; 
                    do
                    {
                        double u1 = rnd();
                        double u2 = rnd();
                        double phi = 2*M_PI*u2;
                        double theta = atan(alpha_g*sqrt(u1/(1-u1)));
                        halfv = normalize( u*cos(phi)*sin(theta) + v*sin(phi)*sin(theta) + orienting_normal*cos(theta) );
                        dir = reflect(nextRay.direction, halfv);
                    }while(dot(dir, orienting_normal) < 0);
                    bsdf = F(dot(-1*nextRay.direction, halfv),nowobjcolor) * G(dir,-1*nextRay.direction,orienting_normal,alpha_g) * D_GGX(halfv, orienting_normal,alpha_g) / 4 / std::abs(dot(-1*nextRay.direction, orienting_normal)) / std::abs( dot(dir, orienting_normal) );
                    pdf = D_GGX(halfv, orienting_normal, alpha_g) * std::abs(dot(halfv, orienting_normal)) / 4 / std::abs(dot(halfv, dir));
                    costerm = std::abs(dot(dir, orienting_normal));
                    is_specular = false;
                    /////////NEE
                    int i = accels.select_light();
                    double pa = 1/accels.light_area;
                    Vec3 x_l = accels.light[i]->randompoint();
                    Ray shadowRay = Ray(hit.hitPos, normalize(x_l - hit.hitPos));
                    double length_xl_x = (x_l-hit.hitPos).length();
                    Hit shadow_hit;
                    if(accels.intersect(shadowRay, shadow_hit, 0, nodes))
                    {
                        if( std::abs(shadow_hit.t - length_xl_x) < 1e-8 && dot(shadowRay.direction,shadow_hit.hitNormal) < 0 && dot(orienting_normal, shadowRay.direction) > 1e-8)
                        {
                            
                            Vec3 nee_halfv = normalize(-1*nextRay.direction + shadowRay.direction);
                            Vec3 nee_dir = shadowRay.direction;

                            double bsdf_pdf = D_GGX(nee_halfv, orienting_normal, alpha_g) * std::abs(dot(nee_halfv, orienting_normal)) / 4 / std::abs(dot(nee_halfv, nee_dir)) * dot(shadowRay.direction, orienting_normal)/length_xl_x/length_xl_x;
                            double mis_weight = pa / (bsdf_pdf + pa);

                            Vec3 _bsdf = F(dot(-1*nextRay.direction, nee_halfv), nowobjcolor) * G(nee_dir, -1*nextRay.direction, orienting_normal, alpha_g) * D_GGX(nee_halfv, orienting_normal, alpha_g) / dot(nee_dir, orienting_normal) / dot(-1*nextRay.direction,orienting_normal) / 4;
                            double G_nee = std::abs(dot(shadow_hit.hitNormal,normalize(hit.hitPos - x_l)))*std::abs(dot(orienting_normal,normalize(x_l - hit.hitPos))) / pow(length_xl_x,2);
                            I = I + mis_weight * weight * accels.light[i]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * _bsdf * G_nee *(1/pa);
                        }
                    }
                    /////////
                    break;
                }
            }
            //////////
            nextRay = Ray(hit.hitPos, dir);//次のRay
            weight = weight * bsdf * costerm * (1 / pdf);//ウェイトの更新 
            if(isNan(weight)) std::cout << "hahaha" << std::endl;
            Prr *= 0.96;//ロシアンルーレットの確率を決める.
        
            if(rnd() >= Prr)//ロシアンルーレットの開始だ!!
            {
                break;
            }
            weight = (1/Prr)*weight;//ウェイトにロシアンルーレットの確率を反映させる.
            if(weight == Vec3()) 
            {
                break;
            }
        }
        else//お空にRayが当たった!
        {
            //I = I + weight * ibl.getColor2(nextRay); //環境マップ
            break;
        }
    }
    return I;
}

#endif