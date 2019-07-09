#ifndef ACCEL_H
#define ACCEL_H
#include <cassert>
#include <deque>
#include <memory>
#include "myrandom.h"
#include "figure.h"
#include "hit.h"
#include "ray.h"
#include "bvh.h"

class Accel
{
  public:
    std::deque<std::shared_ptr<Figure>> shapes;
    std::deque<std::shared_ptr<Figure>> light;
    std::vector<double> light_cdf;
    double light_area = 0;

    Accel(){};

    void add(std::shared_ptr<Figure> p)
    {
        shapes.push_back(p);
        if(p->light != nullptr)
        {
            light.push_back(p);
            light_area += p->area();
        }
    };

    void build_light(){
        assert(light_area > 0);
        int num_lights = light.size();

        light_cdf.resize(num_lights);

        light_cdf.at(0) = light[0]->area()/light_area;
        for(int i = 1; i < num_lights; i++) {
            light_cdf.at(i) =light_cdf[i-1] + light[i]->area()/light_area;
        }
    }

    int select_light() {
        assert(light_area > 0);
        double u = rnd();
        return std::lower_bound(light_cdf.begin(), light_cdf.end(), u) - light_cdf.begin();
    }

    //fixme: わざわざfigureへのポインタを返すのはメモリの無駄使い
    bool intersect(const Ray& ray, Hit& hit, int index, BVHnode* node) {
        bool is_hit_aabb = intersectAABBvsRay(node[index].bbox, ray);

        if(!is_hit_aabb) { //AABBに当たらなかった
            return false;
        }
        else {
            if(node[index].children[1] == -1) { //葉ノード
                Hit hit_each;
                bool is_hit = false;
                for(auto& p : node[index].polygons) {
                    if(p->intersect(ray, hit_each)) {
                        if(hit_each.t < hit.t) {
                            hit = hit_each;
                        }
                        is_hit = true;
                    }
                }
                return is_hit;
            }
            else { //中間ノード

                bool is_hit1 = intersect(ray, hit, node[index].children[0], node);
                bool is_hit2 = intersect(ray, hit, node[index].children[1], node);

                return is_hit2 || is_hit1;
            }
        }
    }
};
#endif
