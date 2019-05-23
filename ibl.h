#ifndef IBL_H
#define IBL_H
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <string>

struct IBL
{
    int width;
    int height;
    float *ibl_data;
    IBL(const std::string &filename)
    {
        int n;
        ibl_data = stbi_loadf(filename.c_str(), &width, &height, &n, 0);
    };
    ~IBL()
    {
        stbi_image_free(ibl_data);
    };

    Vec3 getColor2(const Ray &ray) const
    {

        double phi = std::atan2(ray.direction.z, ray.direction.x) - M_PI;
        double theta = std::acos(ray.direction.y);
        phi = phi - 1.3*M_PI;
        if (phi < 0)
            phi += 2 * M_PI;
        if (phi > 2 * M_PI)
            phi -= 2 * M_PI;

        double u = phi / (2 * M_PI);
        double v = theta / M_PI;

        int i = (int)((1 - u) * width);
        int j = (int)(v * height);
        int adr = 3 * i + 3 * width * j;
        return Vec3(ibl_data[adr], ibl_data[adr + 1], ibl_data[adr + 2]);
    }
};

#endif