#ifndef PSC_AABB_H
#define PSC_AABB_H

#include "vec_math.h"

class cubvh
{
  public:
    int left;
    int right;
    float3 n;
    float3 m;
    float3 conec_up;
    float cos_coner_up;
    float3 conec_down;
    float cos_coner_down;
    int ni;
    int idx[1];
};

//class cuedge
//{
//  public:
//    int v0, t0;
//    int v1, t1;
//};

class cuedge
{
  public:
    int v0, v1;
    int t0, t1;
};




class curay
{
  public:
    float3 u0;
    float3 du;
};

class cups
{
  public:
    float3 origin;
    float3 ev0;
    float3 ev1;
    float3 vvn;

    int eidx;
    int tri_idx;
};

class CUParam
{
  public:
    unsigned char *cbvh;
    float3 *vmap, *nmap, *vnmap;
    cuedge *emap;
    int *vimap, *eimap;
};

class PointI
{
  public:
    int   ei;
    float es;
};


class EdgeJ
{
  public:
    float es0;
    int u0_facing;
    float es1;
    int   ei1;
};

class CUType0
{
  public:
    float3 pt;
    int ti;
};


#endif


