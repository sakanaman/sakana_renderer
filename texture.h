#ifndef TEXTURE_H
#define TEXTURE_H
#include <cmath>
#include "vec3.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "perlin.h"

class Texture
{
  public:
    virtual Vec3 getColor1(double u, double v, const Vec3 &p) const = 0;
};

class UniformTexture : public Texture
{
  public:
    Vec3 col;
    UniformTexture(const Vec3 &col_) : col(col_){};
    virtual Vec3 getColor1(double u, double v, const Vec3 &p) const
    {
        return col;
    };
};

class CheckerTexture : public Texture
{
  public:
    Vec3 color1;
    Vec3 color2;
    CheckerTexture(const Vec3 &color1, const Vec3 &color2) : color1(color1), color2(color2){};

    virtual Vec3 getColor1(double u, double v, const Vec3 &p) const
    {
        double sines = std::sin(10 * p.x) * std::sin(10 * p.y) * std::sin(10 * p.z);
        if (sines < 0)
        {
            return color1;
        }
        else
        {
            return color2;
        }
    }
};

class NoiseTexture : public Texture
{
  public:
    double scale;
    perlin noise;
    NoiseTexture(){};
    NoiseTexture(double sc) : scale(sc){};
    virtual Vec3 getColor1(double u, double v, const Vec3 &p) const
    {
        //return Vec3(1,1,1)*0.5*(1+noise.turb(scale*p));
        //return Vec3(1,1,1)*noise.turb(scale*p);
        return Vec3(1, 1, 1) * 0.5 * (1 + std::sin(scale * p.z + 10 * noise.turb(p)));
    };
};

class ImageTexture : public Texture
{
  public:
    int width;
    int height;
    unsigned char *image_data;
    ImageTexture(const std::string &filename)
    {
        int n;
        image_data = stbi_load(filename.c_str(), &width, &height, &n, 3);
    };
    ~ImageTexture()
    {
        stbi_image_free(image_data);
    };

    virtual Vec3 getColor1(double u, double v, const Vec3 &p) const
    {
        int w = (int)(u * width);
        int h = (int)(v * height);
        int adr = 3 * w + 3 * width * h;
        return Vec3(image_data[adr] / 255.0, image_data[adr + 1] / 255.0, image_data[adr + 2] / 255.0);
    };
};

#endif