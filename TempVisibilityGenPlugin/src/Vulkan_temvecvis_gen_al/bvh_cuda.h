#ifndef BVH_CUDA_H
#define BVH_CUDA_H

#include "aabb.h"
#include "psc_aabb.cuh"

int bvh_size(my_node* root, int addr);
int bvh_memcpy(my_node* root, unsigned char* bvh0, int addr);
void bvh_print(unsigned char* bvh0, int addr, int& m);

extern void cu_ps_bvh_flattened_cone_speed( const CUParam &param, cups *ps,  float3 *CURecord_vis, int bn, int tn);
extern void cu_rt_bvh( const CUParam &param, curay *rays, CUType0 *results, int ns );
extern void cu_initialize_ps(
  int *front_back_edge_buf, 
  int *front_back_edge_offset, 
  cups *&ps_buf,
  int *&front_back_edge_tightpacked, 
  int &front_back_edge_tightpacked_size, 
  int &_front_back_edge_tightpacked_size,
  int *front_back_edge_vinfo,
  const CUParam &param, int n_edge, 
  int n_vertex, int n_update );

void cu_ps_bvh_flattened_cone_speed_g0(
  PointI *&vis_tightpacked, int &vis_tightpacked_size, int &_vis_tightpacked_size, int *vis_tightpacked_vinfo, 
  PointI *vis_buf, int *vis_offset,
  const CUParam &param, cups *ps_buf, int *front_back_edge_vinfo, int nps,
  int n_vertex );

void cu_ps_bvh_flattened_cone_speed_testing(
  PointI *vis_buf, int *vis_offset,
  const CUParam &param, cups *ps_buf, int *front_back_edge_vinfo, int nps, int n_vertex, int n_update);

void my_check_cuda();

#endif 

