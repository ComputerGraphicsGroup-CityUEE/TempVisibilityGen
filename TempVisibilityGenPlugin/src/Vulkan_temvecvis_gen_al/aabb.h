#ifndef AABB_H
#define AABB_H

#include <stdio.h>
#include "g_vector.h"
#include <vector>

class aabb
{
  public:
    int idx;
    FLOAT3 n;  // min 
    FLOAT3 m;  // max 


    aabb()
    {
      idx = 0;
    }

    void creat_box(FLOAT3 vmax, FLOAT3 vmin, int _idx)
    {
      m = vmax;
      n = vmin;
      idx = _idx;
    }

    void print_box()
    {
      printf("max: %f %f %f\n", m.x, m.y, m.z);
      printf("min: %f %f %f\n", n.x, n.y, n.z);
    }
};


class my_node
{
  public:
    my_node *left, *right;
    aabb box;
    std::vector<int> idx;

    FLOAT3 conec_up, conec_down;
    float  coner_up, coner_down;
    float  cos_coner_up, cos_coner_down;

    my_node()
    {
      left = NULL;
      right = NULL;
      coner_up = 0;
      coner_down = 0;
      cos_coner_up = 1;
      cos_coner_down = 1;

      //memset(this, 0, sizeof(*this));
    }
};


#endif


