#include "vec3.h"
#include "image.h"
#include "camera.h"
#include "ray.h"
#include "figure.h"
#include "accel.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <float.h>
#include "myrandom.h"
#include "texture.h"
#include "ibl.h"
#include "perlin.h"
#include "bvh.h"
#include "obj.h"
#include "light.h"
//#define USE_DEBUG
#include "radience.h"

Accel accel;

int main()
{
    IBL ibl("texture_ibl/PaperMill_E_3k.hdr"); //to do: select use or nonuse
    int samples = 10;
    Image img(512,512); //横・縦幅をもつImageのオブジェクトの生成
    Vec3 lookfrom( 0, -1.5, -8.5);//(0,-1.5, -8.5)か(0, -2.5, -6.5)
    //Vec3 lookfrom(0, 8.0 * std::sin(12 * M_PI / 180), -8.0 * std::cos(12 * M_PI / 180));
    Pinhole_Camera cam(lookfrom, normalize(Vec3(0,0,1)));
    //Thin_Lens_Camera cam(lookfrom, normalize(Vec3(0,0,1)), 1.5, Vec3(0,-0.1,-1), 8.0);
    //Thin_Lens_Camera cam(lookfrom, normalize(-1 * lookfrom), 1.5, Vec3(0, 0, 0), 9.0);
    double screen_height = 2.0;
    accel.add(std::make_shared<Tri>(Vec3(-4, -4, 2), Vec3(4, -4, -8), Vec3(-4, -4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //下の地面1
    accel.add(std::make_shared<Tri>(Vec3(-4, -4, 2), Vec3(4, -4, -8), Vec3(4, -4, 2),
                                    std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //下の地面2
    accel.add(std::make_shared<Tri>(Vec3(-4, -4, 2), Vec3(-4, 4, 2), Vec3(-4, -4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.75, 0.25, 0.25)), 0, nullptr)); //左横の面1
    accel.add(std::make_shared<Tri>(Vec3(-4, 4, 2), Vec3(-4, -4, -8), Vec3(-4, 4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.75,0.25,0.25)), 0, nullptr)); //左横の面2
    accel.add(std::make_shared<Tri>(Vec3(4, -4, 2), Vec3(4, 4, 2), Vec3(4, -4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.25, 0.25, 0.75)), 0, nullptr)); //右横の面1
    accel.add(std::make_shared<Tri>(Vec3(4, 4, 2), Vec3(4, -4, -8), Vec3(4, 4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.25,0.25,0.75)), 0, nullptr)); //右横の地面2
    accel.add(std::make_shared<Tri>(Vec3(-4, 4, 2), Vec3(4, 4, -8), Vec3(4, 4, 2),
                                    std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //上の面1
    accel.add(std::make_shared<Tri>(Vec3(-4, 4, 2), Vec3(4, 4, -8), Vec3(-4, 4, -8),
                                    std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //上の面2
     accel.add(std::make_shared<Tri>(Vec3(-4, -4, 2), Vec3(-4, 4, 2), Vec3(4, 4, 2),
                                  std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //奥の面1
     accel.add(std::make_shared<Tri>(Vec3(4, 4, 2), Vec3(4, -4, 2), Vec3(-4, -4, 2),
                                     std::make_shared<UniformTexture>(Vec3(0.75,0.75,0.75)), 0, nullptr)); //奥の面2
    accel.add(std::make_shared<Sphere>(Vec3(2.3,-2.3,-1.5),1.7,
                                        std::make_shared<UniformTexture>(Vec3(0.99,0.99,0.99)),2, nullptr));
    accel.add(std::make_shared<Sphere>(Vec3(-2.3,-2.3, 0.3),1.7,
                                        std::make_shared<UniformTexture>(Vec3(0.99,0.99,0.99)),1, nullptr));
    accel.add(std::make_shared<Tri>(Vec3(-4,3.99,2),Vec3(0,3.99,-0.5),Vec3(-4,3.99,-8),
                                    std::make_shared<UniformTexture>(Vec3()), 0,
                                    std::make_shared<Isotropic_Light>(2)));
    accel.add(std::make_shared<Tri>(Vec3(-4,3.99,-8),Vec3(0,3.99,-0.5),Vec3(4,3.99,-8),
                                    std::make_shared<UniformTexture>(Vec3()), 0,
                                    std::make_shared<Isotropic_Light>(2)));
    accel.add(std::make_shared<Tri>(Vec3(4,3.99,-8),Vec3(0,3.99,-0.5),Vec3(4,3.99,2),
                                    std::make_shared<UniformTexture>(Vec3()), 0,
                                    std::make_shared<Isotropic_Light>(2)));
    accel.add(std::make_shared<Tri>(Vec3(4,3.99,2),Vec3(0,3.99,-0.5),Vec3(-4,3.99,2),
                                    std::make_shared<UniformTexture>(Vec3()), 0,
                                    std::make_shared<Isotropic_Light>(2)));
    // accel.add(std::make_shared<Sphere>(Vec3(0, 4.00, -1.5), 1.5,
    //                                     std::make_shared<UniformTexture>(Vec3()),0,
    //                                     std::make_shared<Isotropic_Light>(5)));
    
    
    //bunnyちゃん
    // TriangleMesh mesh;
    // double scale = 1.0;
    // input_obj("object/bunny.obj", mesh, scale);
    // for(int i = 0; i < mesh.num_shapes; i++)
    // {
    //     Vec3 a = mesh.as[i]-Vec3(-0.2,3.2,1);
    //     Vec3 b = mesh.bs[i]-Vec3(-0.2,3.2,1);
    //     Vec3 c = mesh.cs[i]-Vec3(-0.2,3.2,1);
    //     accel.add(std::make_shared<Tri>(a, b, c,
    //                 std::make_shared<UniformTexture>(Vec3( 0.99, 0.99, 0.99)), 2, nullptr));
    // }
    constructBVH(accel.shapes);
    accel.build_light();
#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
    for (int y = 0; y < img.height; y++) 
    {
        for (int x = 0; x < img.width; x++) 
        {
            const int image_index = y * img.width + x;
            Vec3 accumulated_radiance = Vec3();
            for(int i=0; i < samples; i++)
            {
                double u = screen_height*(2*x - img.width) / (2 * img.height)
                            +screen_height*rnd()/img.height;
                double v = screen_height*(2*y - img.height)/(2*img.height)
                            +screen_height*rnd()/img.height;
                double w;
                Ray ray = cam.getRay(u,v,w);
                //fixme: 引数にiblを直接渡す実装は汚い
                #ifndef USE_DEBUG
                accumulated_radiance = accumulated_radiance + getColor(ray, ibl,accel) / samples;
                #else
                accumulated_radiance = accumulated_radiance + Debug_normal(ray, accel) / samples;
                //accumulated_radiance = accumulated_radiance + Debug_depth(ray, accel) / samples;
                //accumulated_radiance = accumulated_radiance + Debug_ao(ray, accel) / samples;
                #endif
            }
            img.data[image_index] = img.data[image_index] + accumulated_radiance;
        }
        fprintf(stderr, "\r[%3d]", int(100 * y /(img.height - 1))); 
    }
    std::cout << "used node count: "<< used_node_count << endl;
    img.gamma_correction_clamp();
    img.ppm_output();
}