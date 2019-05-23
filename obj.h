#ifndef OBJ_H
#define OBJ_H
#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_DOUBLE
#include "tiny_obj_loader.h"
#include "vec3.h"
#include <iostream>
#include <vector>
#include <string>

struct TriangleMesh
{
    size_t num_shapes;
    std::vector<Vec3> as;
    std::vector<Vec3> bs;
    std::vector<Vec3> cs;
};


void input_obj(const std::string& filename, TriangleMesh &mesh, double scale)
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str(),__null, true, true);
  
    if (!err.empty()) 
    {
        std::cerr << err << std::endl;
    }

    if (!ret)
    {
        exit(1);
    }

    for (size_t s = 0; s < shapes.size(); s++)
    {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            int fv = shapes[s].mesh.num_face_vertices[f];
            for (size_t v = 0; v < fv; v++)
            {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
                tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
                tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
                if(v==0)mesh.as.push_back(scale * Vec3(vx,vy,-vz));
                if(v==1)mesh.bs.push_back(scale * Vec3(vx,vy,-vz));
                if(v==2)mesh.cs.push_back(scale * Vec3(vx,vy,-vz));
            }
            index_offset += fv;
        }
    }
    mesh.num_shapes = mesh.as.size();
}

#endif