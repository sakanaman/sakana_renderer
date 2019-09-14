#ifndef VOLUME_H
#define VOLUME_H
#include "vec3.h"
#include "myrandom.h"

double sample_scattering_point(const bool homo, const double s_max, const double sigma_e){
    if(homo){
        double u = rnd();
        double s = -1.0 * std::log(1 - u) / sigma_e;
        return s;
    }
    else {
        //TODO : implement ray marching
        return 10000000;
    }
}

Vec3 phase_sampling(int scatter_type, double* pdf) {
    if(scatter_type == 0){ //equal party
        Vec3 result = randomsphere();
        *pdf= 1.0/4/M_PI;
        return result;
    }
    return Vec3(); //TODO: implement several type of scattering
}

double phase_funtion(int scatter_type) {
    if(scatter_type == 0) {
        return 1.0/4/M_PI;
    }
    return 1.0;  //TODO: implement several type of scattering
}

double Transmit(bool homo, double sigma_e, double l) {
    if(homo) {
        return std::exp(-1.0 * sigma_e * l);
    }
    return 1.0;//todo
}
#endif