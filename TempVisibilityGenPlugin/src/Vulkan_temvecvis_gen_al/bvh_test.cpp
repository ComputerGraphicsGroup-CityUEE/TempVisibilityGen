#include <bvh/bvh.hpp>
#include <bvh/vector.hpp>
#include <bvh/triangle.hpp>
#include <bvh/ray.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/locally_ordered_clustering_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>

#include "bvh_test.h"


using Scalar   = float;
using Vector3  = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Ray      = bvh::Ray<Scalar>;
using Bvh      = bvh::Bvh<Scalar>;
using Morton = uint32_t;

Bvh bvh0;
vector<Triangle> bvh0_triangles;
bvh::LocallyOrderedClusteringBuilder<Bvh, Morton> builder(bvh0);
cubvh *cubvh_buffer = NULL;
int cubvh_buffer_size = 0;

#include "glm/glm.hpp"
//#include"obj2gltf.h"


int bvh0_node_count()
{
  return bvh0.node_count;
}

#include"g_pfm.h"
void save_cubvhX(cubvh* bvhbuf, int nb)
{
  int w = 512;
  int h = nb/w +1;

  printf("bvh, nb, w, h %d %d %d\n",  nb, w, h );
  GPfm intbuf, maxbuf, minbuf;

  intbuf.load(w, h);
  maxbuf.load(w, h);
  minbuf.load(w, h);

  for(int i=0; i<w*h; i++)
  {
    if(i<nb)
    {
      intbuf.fm[i] =FLOAT3(bvhbuf[i].left, bvhbuf[i].right, bvhbuf[i].ni);
      maxbuf.fm[i] =FLOAT3(bvhbuf[i].m.x, bvhbuf[i].m.y, bvhbuf[i].m.z);
      minbuf.fm[i] =FLOAT3(bvhbuf[i].n.x, bvhbuf[i].n.y, bvhbuf[i].n.z);
    }
  }
  intbuf.save("intbuf.pfm");
  maxbuf.save("maxbuf.pfm");
  minbuf.save("minbuf.pfm");
}


cubvh* bvh0_to_cubvh()
{
  if( cubvh_buffer_size < bvh0.node_count )
  {
    if(cubvh_buffer)
      free(cubvh_buffer);
    cubvh_buffer_size = int( bvh0.node_count * 1.5 );
    cubvh_buffer = (cubvh*) malloc( cubvh_buffer_size*sizeof(cubvh) );
  }
  int i;

  memset(cubvh_buffer, 0, bvh0.node_count*sizeof(cubvh) );

  //int w = 512;
  //int h = bvh0.node_count/w+1;
  //printf("bvh0.node_count, w, h, %d %d %d\n",bvh0.node_count, w, h);
  //GPfm ii, lr, mm, nn;
  //ii.load(w, h);
  //lr.load(w, h);
  //mm.load(w, h);
  //nn.load(w, h);

  for( i=0; i<bvh0.node_count; i++ )
  {
    cubvh &t = cubvh_buffer[i];
    t.m = make_float3(bvh0.nodes[i].bounds[1], bvh0.nodes[i].bounds[3], bvh0.nodes[i].bounds[5]);
    t.n = make_float3(bvh0.nodes[i].bounds[0], bvh0.nodes[i].bounds[2], bvh0.nodes[i].bounds[4]);

  //mm.fm[i] = FLOAT3(bvh0.nodes[i].bounds[1], bvh0.nodes[i].bounds[3], bvh0.nodes[i].bounds[5]);
  //nn.fm[i] = FLOAT3(bvh0.nodes[i].bounds[0], bvh0.nodes[i].bounds[2], bvh0.nodes[i].bounds[4]);

    if( bvh0.nodes[i].is_leaf )
    {
      t.ni = 1;
      t.idx[0] = int( bvh0.primitive_indices[bvh0.nodes[i].first_child_or_primitive] );

      //ii.fm[i] = FLOAT3(t.ni, t.idx[0], 1.0);
      //lr.fm[i] = FLOAT3(0,0, 1.0);
    }
    else
    {
      t.left = int( (bvh0.nodes[i].first_child_or_primitive + 0) * sizeof(cubvh) );
      t.right = int( (bvh0.nodes[i].first_child_or_primitive + 1) * sizeof(cubvh) );

      //ii.fm[i] = FLOAT3(0, 0, 1.0);
      //lr.fm[i] = FLOAT3(t.left,t.right, 1.0);

    }
  }
  //save_cubvhX(cubvh_buffer, bvh0.node_count);


  //ii.save("ii.pfm");
  //lr.save("lr.pfm");
  //mm.save("mm.pfm");
  //nn.save("nn.pfm");
  return cubvh_buffer;
}

void bvh0_to_cubvhVec(std::vector<cubvh> &dat)
{
  printf("DLL 00 cubvh_buffer_size %d\n ", cubvh_buffer_size);
  printf("DLL 00 dat.size() %d\n ", dat.size());
 if(cubvh_buffer_size < bvh0.node_count)
  {
    cubvh_buffer_size = int( bvh0.node_count * 1.5 );
    dat.clear();
    dat.resize(cubvh_buffer_size);
    printf("DLL 01 cubvh_buffer_size %d\n ", cubvh_buffer_size);
    printf("DLL 01 dat.size() %d\n ", dat.size());
  }
    
  int i;
  for(i=0; i<bvh0.node_count; i++)
  {
    cubvh t;
    t.m = make_float3(bvh0.nodes[i].bounds[1], bvh0.nodes[i].bounds[3], bvh0.nodes[i].bounds[5]);
    t.n = make_float3(bvh0.nodes[i].bounds[0], bvh0.nodes[i].bounds[2], bvh0.nodes[i].bounds[4]);
    if( bvh0.nodes[i].is_leaf )
    {
      t.ni = 1;
      t.idx[0] = int( bvh0.primitive_indices[bvh0.nodes[i].first_child_or_primitive] );
    }
    else
    {
      t.left = int( (bvh0.nodes[i].first_child_or_primitive + 0) * sizeof(cubvh) );
      t.right = int( (bvh0.nodes[i].first_child_or_primitive + 1) * sizeof(cubvh) );
    }
    dat[i] = t;
  }

}


void bvh0_update( const MyVertex *vbuf, const std::vector<uint32_t> &ibuf )
{

  int n_face = int(ibuf.size()/3);
  int i;

  //GPfm v0m, v1m, v2m;
  //GPfm vn0m, vn1m, vn2m;
  //int w = 512;
  //int h = n_face/w +1;

  //v0m.load(w, h);
  //v1m.load(w, h);
  //v2m.load(w, h);
  //vn0m.load(w, h);
  //vn1m.load(w, h);
  //vn2m.load(w, h);

  if(bvh0_triangles.size()!=n_face)
    bvh0_triangles = vector<Triangle>(n_face);
  //#pragma omp parallel for
  for( i=0; i<n_face; i++ )
  {
    Vector3 &v0 = *((Vector3*)&vbuf[ibuf[3*i+0]]);
    Vector3 &v1 = *((Vector3*)&vbuf[ibuf[3*i+1]]);
    Vector3 &v2 = *((Vector3*)&vbuf[ibuf[3*i+2]]);
    bvh0_triangles[i] = Triangle(v0,v1,v2);

    //v0m.fm[i] = FLOAT3(v0[0], v0[1], v0[2]);
    //v1m.fm[i] = FLOAT3(v1[0], v1[1], v1[2]);
    //v2m.fm[i] = FLOAT3(v2[0], v2[1], v2[2]);
    //vn0m.fm[i] = FLOAT3(v0[3], v0[4], v0[5]);
    //vn1m.fm[i] = FLOAT3(v1[3], v1[4], v1[5]);
    //vn2m.fm[i] = FLOAT3(v2[3], v2[4], v2[5]);


    //if(i%100==0)
    //{
    //  printf("i,\tibuf[3*i+0], %d\t%d\t\n", i, ibuf[3*i+0]);
    //  printf("v0\t%f\t%f\t%f\n",v0[0], v0[1], v0[2]);
    //  printf("vbuf[ibuf[3*i+0]]\t%f\t%f\t%f\n",vbuf[ibuf[3*i+0]].pos.x,vbuf[ibuf[3*i+0]].pos.z,vbuf[ibuf[3*i+0]].pos.z);
    //}
  }
  //v0m.save("v0m.pfm");  
  //v1m.save("v1m.pfm");  
  //v2m.save("v2m.pfm");  
  //vn0m.save("vn0m.pfm");  
  //vn1m.save("vn1m.pfm");  
  //vn2m.save("vn2m.pfm");  

  auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(bvh0_triangles.data(), bvh0_triangles.size());
  auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), bvh0_triangles.size());
  builder.build(global_bbox, bboxes.get(), centers.get(), bvh0_triangles.size());
}

CloseHitResult* bvh0_raytrace(CloseHitResult *res, FLOAT3 ray_origin, FLOAT3 ray_dir)
{
  Ray ray(
    Vector3(ray_origin.x, ray_origin.y, ray_origin.z), // origin
    Vector3(ray_dir.x, ray_dir.y, ray_dir.z), // direction
    0.0,                    // minimum distance
    100.0                   // maximum distance
  );
  bvh::ClosestPrimitiveIntersector<Bvh, Triangle> primitive_intersector(bvh0, bvh0_triangles.data());
  bvh::SingleRayTraverser<Bvh> traverser(bvh0);

  std::optional<bvh::ClosestPrimitiveIntersector<Bvh, Triangle>::Result> hit;
  hit = traverser.traverse(ray, primitive_intersector);
  if(hit)
  {
    res->tidx = hit->primitive_index;
    res->bc = FLOAT3( 1-hit->intersection.u-hit->intersection.v, hit->intersection.u, hit->intersection.v );
    return res;
  }else
    return NULL;

  //if (hit) 
  //{
  //    auto triangle_index = hit->primitive_index;
  //    auto intersection = hit->intersection;
  //    std::cout << "Hit triangle " << triangle_index << "\n"
  //              << "distance: "    << intersection.t << "\n"
  //              << "u: "           << intersection.u << "\n"
  //              << "v: "           << intersection.v << std::endl;
  //}
}

