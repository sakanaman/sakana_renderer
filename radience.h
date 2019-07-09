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

class PathInfo {
public:
    double pdf_w;
    double costerm;
    Vec3 f;
    Vec3 dir;
    Vec3 throughput;
    Vec3 I;
    bool is_specular;
    PathInfo() {
        pdf_w = 1.0;
        costerm = 1.0;
        f = dir = Vec3(0,0,0);
        is_specular = false;
        throughput = Vec3(1.0, 1.0, 1.0);
        I = Vec3(0,0,0);
    }
};

//
//~Process of Path Tracing~
//while(ray unterminate) {
// caluculate color with last ray
// russian roulette
// if(ok_next) 
// {
//   bsdf sampling
//   light sampling
// {
// else
// { 
//   terminate
// }
//}
Vec3 getColor(const Ray& ray, const IBL& ibl, Accel& accel) 
{
    Ray nextRay = ray;
    PathInfo path;
    double Prr = 1.0;
    for(int depth = 0;; depth++) 
    {
        Hit hit;
        bool is_hit = accel.intersect(nextRay, hit, 0, nodes);
        if(is_hit) 
        {
            double cosine = dot(hit.hitNormal, nextRay.direction);
            if(hit.hitShape->light != nullptr && cosine < 0) // hit light surface
            { 
                cosine = std::abs(cosine);

                if(path.is_specular == true) // traced specular surface
                { 
                    path.I = path.I + path.throughput * hit.hitShape->light->radi_intencity(cosine);
                    break; // terminate ray after hitting light
                }
                if(depth == 0) { // when path length is 0
                    path.I = hit.hitShape->light->radi_intencity(cosine);
                    break; // terminate ray after hitting light
                }

                // evaluate MIS contribution (from BRDF sampling)
                double l = (hit.hitPos - nextRay.origin).length();
                double light_pdf = 1/accel.light_area * l * l / path.costerm;
                double mis_weight = path.pdf_w/(light_pdf + path.pdf_w);
                
                path.I = path.I + mis_weight * path.throughput * hit.hitShape->light->radi_intencity(cosine);
                break;
            }

            //russian roulette
            double Prr = std::max(path.throughput.x, std::max(path.throughput.y, path.throughput.z));
            if(path.is_specular) {
                Prr = 0.96;
            }
            if(rnd() >= Prr) {
                break;
            }
            path.throughput = path.throughput * 1/Prr;

            ///BSDF Sampling & Light Sampling
            const Vec3 orienting_normal = dot(-1*nextRay.direction,hit.hitNormal)>0 ? hit.hitNormal : -1.0*hit.hitNormal;
            const Vec3 nowobjcolor = hit.hitShape->texture->getColor1(hit.u,hit.v,hit.hitPos);
            switch(hit.hitShape->material) {
                case 0: //diffuse
                {
                    //BRDF sampling
                    path.f = nowobjcolor * (1/M_PI);
                    path.dir = randomCosineHemisphere(orienting_normal);
                    path.costerm = std::abs(dot(path.dir,orienting_normal));
                    path.pdf_w = path.costerm/M_PI;
                    path.is_specular = false;

                    //Light Sampling
                    int light_id = accel.select_light();
                    double pdf_area = 1/accel.light_area;
                    Vec3 x_l = accel.light[light_id]->randompoint();
                    double length_xl_x = (x_l-hit.hitPos).length();

                    Ray shadowRay(hit.hitPos, normalize(x_l - hit.hitPos));
                    Hit shadow_hit;
                    if(accel.intersect(shadowRay, shadow_hit, 0, nodes)) 
                    {
                        if(std::abs(shadow_hit.t - length_xl_x) < 1e-8 && dot(shadowRay.direction,shadow_hit.hitNormal) < 0 && dot(orienting_normal, shadowRay.direction) > 1e-8)
                        {
                            double bsdf_pdf = dot(shadowRay.direction, orienting_normal)/M_PI * dot(shadowRay.direction, orienting_normal)/length_xl_x/length_xl_x;
                            double mis_weight = pdf_area/(bsdf_pdf + pdf_area);

                            double G = std::abs(dot(shadow_hit.hitNormal,normalize(hit.hitPos - x_l)))*std::abs(dot(orienting_normal,normalize(x_l - hit.hitPos))) / pow(length_xl_x,2);
                            //I = I +  mis_weight * weight * accels.light[i]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * bsdf * G *(1/pa);
                            path.I = path.I +  mis_weight * path.throughput * accel.light[light_id]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * path.f * G *(1/pdf_area);
                        }
                    }
                    break;
                }
                case 1: // mirror
                {
                    path.dir = reflect(nextRay.direction,orienting_normal);
                    path.costerm = std::abs(dot(path.dir,orienting_normal));
                    path.f = nowobjcolor*(1/path.costerm);
                    path.pdf_w = 1.0;
                    path.is_specular = true;
                    break;
                }
                case 2: // dielectric
                {
                    path.is_specular = true;
                    double u = rnd();
                    bool in = (dot(hit.hitNormal, orienting_normal) > 0);
                    double nt = in? 1.5 : 1.0;
                    double ni = in? 1.0 : 1.5;
                    double costhetaI = dot(orienting_normal, -1 * nextRay.direction);
                    double sinthetaT = ni/nt * std::sqrt(std::max(1 - costhetaI * costhetaI,0.0));
                    if(sinthetaT >= 1) { // total reflect
                        path.pdf_w = 1.0;
                        path.dir = reflect(nextRay.direction, orienting_normal);
                        //if(std::abs( dot(dir, orienting_normal) - dot(orienting_normal, -1*nextRay.direction)) >= 1e-8) std::cout << "fail reflect" << std::endl;
                        path.costerm = std::abs(dot(path.dir, orienting_normal));
                        path.f = nowobjcolor / path.costerm;
                    }
                    else {
                        double costhetaT = std::sqrt(1 - sinthetaT*sinthetaT);
                        double Fr = F_die(costhetaI, 1.0 , 1.5);
                        double u = rnd();
                        if(u < Fr) {
                            path.pdf_w = Fr;
                            path.dir = reflect(nextRay.direction, orienting_normal);
                            path.costerm = std::abs(dot(path.dir, orienting_normal));
                            path.f = Fr * nowobjcolor / path.costerm;
                        }
                        else {
                            path.pdf_w = 1.0 - Fr;
                            path.dir = refract(-1*nextRay.direction, orienting_normal, costhetaI, costhetaT, ni/nt);
                            path.costerm = std::abs(dot(path.dir, orienting_normal));
                            path.f = (ni*ni)/(nt*nt) * (1.0 - Fr) / path.costerm * nowobjcolor;
                        }
                    }
                    break;
                }
                case 3: //GGX reflection (walter)
                {
                    ///BRDF Sampling
                    double alpha_g = 0.6 * 0.6;
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
                        path.dir = reflect(nextRay.direction, halfv);
                    }while(dot(path.dir, orienting_normal) < 0);

                    path.f = F(dot(-1*nextRay.direction, halfv),nowobjcolor) * G(path.dir,-1*nextRay.direction,orienting_normal,alpha_g) * D_GGX(halfv, orienting_normal,alpha_g) / 4 / std::abs(dot(-1*nextRay.direction, orienting_normal)) / std::abs( dot(path.dir, orienting_normal) );
                    path.pdf_w = D_GGX(halfv, orienting_normal, alpha_g) * std::abs(dot(halfv, orienting_normal)) / 4 / std::abs(dot(halfv, path.dir));
                    path.costerm = std::abs(dot(path.dir, orienting_normal));
                    path.is_specular = false;

                    /////////Light Sampling
                    int i = accel.select_light();
                    double pa = 1/accel.light_area;
                    Vec3 x_l = accel.light[i]->randompoint();

                    Ray shadowRay = Ray(hit.hitPos, normalize(x_l - hit.hitPos));
                    double length_xl_x = (x_l-hit.hitPos).length();
                    Hit shadow_hit;

                    if(accel.intersect(shadowRay, shadow_hit, 0, nodes))
                    {
                        if( std::abs(shadow_hit.t - length_xl_x) < 1e-8 && dot(shadowRay.direction,shadow_hit.hitNormal) < 0 && dot(orienting_normal, shadowRay.direction) > 1e-8)
                        {
                            
                            Vec3 nee_halfv = normalize(-1*nextRay.direction + shadowRay.direction);
                            Vec3 nee_dir = shadowRay.direction;

                            double bsdf_pdf = D_GGX(nee_halfv, orienting_normal, alpha_g) * std::abs(dot(nee_halfv, orienting_normal)) / 4 / std::abs(dot(nee_halfv, nee_dir)) * dot(shadowRay.direction, orienting_normal)/length_xl_x/length_xl_x;
                            double mis_weight = pa / (bsdf_pdf + pa);

                            Vec3 _bsdf = F(dot(-1*nextRay.direction, nee_halfv), nowobjcolor) * G(nee_dir, -1*nextRay.direction, orienting_normal, alpha_g) * D_GGX(nee_halfv, orienting_normal, alpha_g) / dot(nee_dir, orienting_normal) / dot(-1*nextRay.direction,orienting_normal) / 4;
                            double G_nee = std::abs(dot(shadow_hit.hitNormal,normalize(hit.hitPos - x_l)))*std::abs(dot(orienting_normal,normalize(x_l - hit.hitPos))) / pow(length_xl_x,2);
                            // I = I + mis_weight * weight * accels.light[i]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * _bsdf * G_nee *(1/pa);
                            path.I = path.I + mis_weight * path.throughput * accel.light[i]->light->radi_intencity(-1*dot(shadowRay.direction,shadow_hit.hitNormal)) * _bsdf * G_nee *(1/pa);
                        }
                    }
                    
                    break;
                }  
            }
            nextRay  = Ray(hit.hitPos, path.dir);

            path.throughput = path.throughput * path.f * path.costerm * (1/path.pdf_w);

            if(isNan(path.throughput)) std::cout << "nan include in weight" << std::endl;
            if(path.throughput.x < 1e-8 && path.throughput.y < 1e-8 && path.throughput.z < 1e-8) 
            {
                break;
            }
        }
        else // hit sky
        {
            break;
        }
    }
    return path.I;
}

#endif