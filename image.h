#ifndef IMAGE_H
#define IMAGE_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "vec3.h"

class Image
{
  public:
    int width;
    int height;
    Vec3 *data;

    Image(int _width, int _height)
    {
        width = _width;
        height = _height;
        data = new Vec3[width * height];
    };

    ~Image()
    {
        delete[] data;
    };


    void ppm_output() const
    {
        std::ofstream file("output.ppm");
        file << "P3" << std::endl;
        file << width << " " << height << std::endl;
        file << "255" << std::endl;
        for (int i = 0; i < height * width; i++)
        {
            Vec3 color_screen = 255 * data[i];
            int r = (int)color_screen.x;
            int g = (int)color_screen.y;
            int b = (int)color_screen.z;
            file << r << " " << g << " " << b << std::endl;
        }
        file.close();
    };

    void gamma_correction_clamp()
    {
        for (int i = 0; i < width * height; i++)
        {
            data[i].x = data[i].x <= 0 ? 0 : std::min(std::pow(data[i].x, 1.0 / 2.2), 1.0);
            data[i].y = data[i].y <= 0 ? 0 : std::min(std::pow(data[i].y, 1.0 / 2.2), 1.0);
            data[i].z = data[i].z <= 0 ? 0 : std::min(std::pow(data[i].z, 1.0 / 2.2), 1.0);
        }
    };
};

#endif
