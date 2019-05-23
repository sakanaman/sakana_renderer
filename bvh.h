#ifndef BVH_H
#define BVH_H
#include<deque>
#include <memory>
#include <algorithm>
#include "figure.h"

double INF = 999999;
double T_tri = 1;
double T_aabb = 1;
int used_node_count = 0;

class BVHnode
{
    public:
        double bbox[2][3];
        int children[2];
        std::deque<std::shared_ptr<Figure>> polygons;
};

BVHnode nodes[100000];

double surfaceArea(const double bbox[2][3])
{
    double dx = bbox[1][0] - bbox[0][0];
    double dy = bbox[1][1] - bbox[0][1];
    double dz = bbox[1][2] - bbox[0][2];
    return 2*(dx*dy + dy*dz + dz*dx);
}

void emptyAABB(double bbox[2][3])
{   
    bbox[0][0] = bbox[0][1] = bbox[0][2] = INF;
    bbox[1][0] = bbox[1][1] = bbox[1][2] = -1 * INF;
    
}

void mergeAABB(const double bbox1[2][3], const double bbox2[2][3], double result[2][3])
{
    for (int j = 0; j < 3; j++)
    {
        result[0][j] = std::min(bbox1[0][j], bbox2[0][j]);
        result[1][j] = std::max(bbox1[1][j], bbox2[1][j]);
    }
}

void createAABBfromTriangles(const std::deque<std::shared_ptr<Figure>> &triangles, double bbox[2][3])
{   
    emptyAABB(bbox);
    for(auto& p : triangles) {
        mergeAABB(p->bbox, bbox, bbox);
    }
}

void makeLeaf(std::deque<std::shared_ptr<Figure>> &polygons, BVHnode *node)
{
    node->children[0] = node->children[1] = -1;
    node->polygons = polygons;
}

void constructBVH_internal (std::deque<std::shared_ptr<Figure>> &polygons, int nodeIndex) 
{
    BVHnode* node = &nodes[nodeIndex];
    createAABBfromTriangles(polygons, node->bbox);
    double bestCost = T_tri * polygons.size();
    int bestAxis = -1;
    int bestSplitIndex = -1;
    double SA_root = surfaceArea(node->bbox);
    for (int axis=0; axis<3; axis++)
    {
        sort(polygons.begin(), polygons.end(),
            [axis](const std::shared_ptr<Figure> a, const std::shared_ptr<Figure> b)
            {
                return a->bbox[1][axis] + a->bbox[0][axis] < b->bbox[1][axis] + b->bbox[0][axis];  
            });
        std::deque<std::shared_ptr<Figure>> s1;
        std::deque<std::shared_ptr<Figure>> s2(polygons);
        double s1bbox[2][3];
        emptyAABB(s1bbox);
        std::deque<double> s1SA(polygons.size()+1, INF);
        std::deque<double> s2SA(polygons.size()+1, INF);

        for (int i=0; i<=polygons.size(); i++) {
            s1SA[i] = fabs(surfaceArea(s1bbox));
            if(s2.size() > 0)
            {
                std::shared_ptr<Figure> p = s2.front();
                s1.push_back(p);
                s2.erase(s2.begin());
                mergeAABB(s1bbox, p->bbox, s1bbox);
            }
        }
        double s2bbox[2][3];
        emptyAABB(s2bbox);
        for (int i=polygons.size(); i>=0; i--)
        {
            s2SA[i] = fabs(surfaceArea(s2bbox));
            if (s1.size() > 0 && s2.size() > 0)
            {
                double cost = 2*T_aabb + (s1SA[i]*s1.size() + s2SA[i]*s2.size()) * T_tri / SA_root;

                if (cost < bestCost)
                {
                    bestCost = cost;
                    bestAxis = axis;
                    bestSplitIndex = i;
                }
            }
            if (s1.size() > 0)
            {
                std::shared_ptr<Figure> p = s1.back();
                s2.push_front(p);
                s1.pop_back();
                mergeAABB(s2bbox, p->bbox, s2bbox);
            }
        }
    }
    if(bestAxis ==  -1)
    {
        makeLeaf(polygons, node);
    }
    else
    {
        sort(polygons.begin(), polygons.end(),
            [bestAxis](const std::shared_ptr<Figure> a, const std::shared_ptr<Figure> b){
                return a->bbox[0][bestAxis] + a->bbox[1][bestAxis] < b->bbox[0][bestAxis] + b->bbox[1][bestAxis];
            });
        used_node_count++;
        node->children[0] = used_node_count;
        used_node_count++;
        node->children[1] = used_node_count;

        std::deque<std::shared_ptr<Figure>> left(polygons.begin(), polygons.begin() + bestSplitIndex);
        std::deque<std::shared_ptr<Figure>> right(polygons.begin() + bestSplitIndex, polygons.end());

        constructBVH_internal(left, node->children[0]);
        constructBVH_internal(right, node->children[1]);
    }
}

void constructBVH(std::deque<std::shared_ptr<Figure>> &polygons)
{
    used_node_count = 0;
    constructBVH_internal(polygons, 0);
}



#endif