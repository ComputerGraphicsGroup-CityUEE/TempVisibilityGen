#ifndef CAL_VISIBILITY_CUDA_H
#define CAL_VISIBILITY_CUDA_H

#include <vector>
//#include "my_vkgltfmodel.h"
#include "bvh_test.h"
#include "psc_aabb.cuh"

//struct GEdge
//{
//  int v0, t0;
//  int v1, t1;
//};
#include "shaders/my_raytrace_datatypes.h"

class VisGen
{
  public:
    CUParam vis_param;
    float3 *vis_record;
    int vis_bvh0_node_count, _vis_bvh0_node_count;

    int *front_back_edge_buf;     // known size
    int *front_back_edge_offset;  // known size
    int *front_back_edge_vinfo;   // known size

    int  front_back_edge_tightpacked_size, _front_back_edge_tightpacked_size;
    int *front_back_edge_tightpacked;    // depends on front_back_edge_tightpacked_size
    cups *ps_buf;                        // depends on front_back_edge_tightpacked_size

    int nps, _nps;
    int nps_update;  // the number of updated line sample, JUST FOR REPORTING ONLY!!!
    PointI *vis_buf;
    int *vis_offset;

    int total_line_sample; 

    PointI *vis_tightpacked;
    int vis_tightpacked_size, _vis_tightpacked_size;
    int *vis_tightpacked_vinfo;

    PointI *my_dat_buf; // cpu
    int *my_dat_vi;     // cpu
    int my_dat_size, _my_dat_size;


    int max_front_back_edge_buf;
    int n_vertex_perpass;

    VisGen();
    ~VisGen();

    void init_cuda_memory(
      int n_vertex,
      int *face_vidx, int *face_eidx, int n_face,
      GEdge* s_edge, int n_edge
    );
    void free_cuda_memory();

    void update_cuda_memory(
      MyVertex *vbuf0, const std::vector< std::vector<uint32_t> > &vidx_distinct,
      FLOAT3 *fnbuf, int n_face, int n_edge
    );
    void update_cuda_memoryX(
      MyVertex *vbuf0, const std::vector<int> &vidx_distinct_all,
      FLOAT3 *fnbuf, int n_face, int n_edge
    );
    void update_bvh( cubvh *dat, int node_count );
    void update_bvh(std::vector<cubvh> &dat, int node_count);
    void gen_vis_rendering3( int n_edge, int vertex0, int n_vertex );

    void testing( int n_edge, int n_vertex, int n_update );

  private:
    VisGen( const VisGen &a );
    VisGen& operator=( const VisGen &a );
};



#endif