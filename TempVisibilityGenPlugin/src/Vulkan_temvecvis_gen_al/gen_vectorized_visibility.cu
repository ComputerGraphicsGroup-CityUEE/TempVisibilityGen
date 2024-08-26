#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include "vec_math.h"
#include "psc_aabb.cuh"
#include "psc_math.cuh"


#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include "gen_vectorized_visibility.h"

#define	G_PI 3.14159265358979323846f

__device__ float flt_epslion = 1.192092896e-07F;
__device__ float flt_max = 3.402823466e+38F;


#define f3abs(a) ( make_float3(fabsf(a.x),fabsf(a.y),fabsf(a.z)) )
//template<typename T> 
//__device__  float3 f3abs( T a )
//{
//  return make_float3( fabsf(a.x), fabsf(a.y), fabsf(a.z) );
//}
__device__ float length2(const float3& v)
{
  return dot(v, v);
}


__device__ void get_coordinate_system(const cups &ps, float3 &ax, float3 &ay, float3 &az)
{
  ay = ps.vvn;
  ax = make_float3(1,0,0);

  //if( fabsf(dot(ax,ay))>.999999 )
  //  ax = make_float3(0,1,0);
  //az = normalize( cross(ax,ay) );

  float3 ref = lerp( make_float3(ay.z, 0.f, -ay.x), ax, abs(ay.y) );
  az = normalize( cross(ref,ay) );

  ax = normalize( cross(ay,az) );
}

__host__ void my_check_cuda();
#include "cu_initialize_ps.inl"


__device__ float3 es_to_float3( float es, const float3 &ev0, const float3 &ev1  )
{
  return es*(ev1-ev0)+(ev1+ev0)/2; 
}
//__device__ float3 es_to_float3( const CUParam &param, float es, int ei )
//{
//  float3 ev0 = param.vmap[param.emap[ei].v0];
//  float3 ev1 = param.vmap[param.emap[ei].v1];
//  return es_to_float3(es, ev0, ev1);
//}

__device__ float LsPlXRatio( float3 v0, float3 fn, float3 ev0, float3 ev1 )
{
  return -dot((ev1+ev0)/2-v0,fn)/dot(ev1-ev0,fn);
}

__device__ float LsPsXRatio( const CUParam &param, int ei0, const float3 &ps_origin, const float3 &ps_ev0, const float3 &ps_ev1 )
{
  float3 au, av, fn;
  au = param.vmap[param.emap[ei0].v0];
  av = param.vmap[param.emap[ei0].v1];
  fn = cross(au - ps_origin, av - ps_origin);
  return -dot((ps_ev1+ps_ev0)/2-ps_origin,fn)/dot(ps_ev1-ps_ev0,fn);
}




__device__ PointI extrude_endpointX( float3 e0, float3 e1, float3 ax, float3 ay, float3 az )
{
  float3 un = normalize(e1-e0);
  PointI pp;
  pp.ei = -1;
  pp.es = acosf( clamp( dot(un,az), -1.f, 1.f ));
  if(dot(un,ax)<0.f)
    pp.es = -pp.es;
  return pp;
}
__device__ PointI getHemiPointI_ev0( const cups &ps, float3 ax, float3 ay, float3 az )
{
  return extrude_endpointX(ps.origin, ps.ev0, ax,ay,az );
}
__device__ PointI getHemiPointI_ev1( const cups &ps, float3 ax, float3 ay, float3 az )
{
  return extrude_endpointX(ps.origin, ps.ev1, ax,ay,az );
}
__device__ PointI getHemiPointI_es( const cups &ps, float es0, float3 ax, float3 ay, float3 az
){
  float3 u1 = es0*(ps.ev1-ps.ev0) + (ps.ev0+ps.ev1)/2;
  return extrude_endpointX(ps.origin, u1, ax, ay, az );
}




//__device__ PointI extrude_endpoint( const CUParam &param, float3 e0, PointI pp, float3 ax, float3 ay, float3 az )
//{
//  float3 ev0 = param.vmap[ param.emap[pp.ei].v0 ];
//  float3 ev1 = param.vmap[ param.emap[pp.ei].v1 ];
//  float3 e1 = pp.es*(ev1-ev0) + (ev0 + ev1)/2;
//  return extrude_endpointX(e0, e1, ax, ay, az );
//}
//__device__ PointI extrude_endpoint( const CUParam &param, PointI pp, float3 e1, float3 ax, float3 ay, float3 az )
//{
//  float3 ev0 = param.vmap[ param.emap[pp.ei].v0 ];
//  float3 ev1 = param.vmap[ param.emap[pp.ei].v1 ];
//  float3 e0 = pp.es*(ev1-ev0) + (ev0 + ev1)/2;
//  return extrude_endpoint(e0, e1, ax, ay, az );
//}
//__device__ float3 extrude_endpoint( float3 e0, float3 e1, float pos )
//{
//  return normalize(e1-e0)*pos + e0;
//}





/////////////////////////////////////////////////////////
//

__device__ bool ps_edge_intersectxa(const CUParam &param, const float3 &vvn, const float3 &ps_origin, const float3 &ps_ev0, const float3 &ps_ev1, int ei0, int ei1,
  
  int ei1v2, 
  float3 &pt,
  float &es1
)
{
  const float3 &pev0 = ps_ev0;
  const float3 &pev1 = ps_ev1;
  const float3 &pt_ev0 = param.vmap[param.emap[ei1].v0];
  const float3 &pt_ev1 = param.vmap[param.emap[ei1].v1];

  //return ps_edge_intersect(ps_origin, pev0, pev1, pt_ev0, pt_ev1, out );

  float3 au = pev0 - ps_origin;
  float3 av = pev1 - ps_origin;
  float3 fn = (cross(au,av) );
  float3 e0 = pt_ev0 - ps_origin;
  float3 e1 = pt_ev1 - ps_origin;
  float3 ee = pt_ev1 - pt_ev0;
  if( fabsf(dot(ee,fn))<flt_epslion || dot(e0,fn)>0 == dot(e1,fn)>0 )
    return false;

  es1 = -dot((pt_ev1+pt_ev0)/2-ps_origin,fn)/(dot(ee,fn));
  es1 = clamp(es1,-.5f,.5f);
  pt = es1*ee + (pt_ev1+pt_ev0)/2;
  // float3 pt = -dot(e0,fn)/dot(ee,fn)*ee + pt_ev0;
  float3 qt = pt - ps_origin;

  if(dot(cross(au,qt),fn)>0 && dot(cross(qt,av),fn)>0 )
  {
    return true;
  }else
    return false;
}



__device__ bool ps_aabb_cone_intersect(const cubvh *ccc, float3 vertex0, float3 vertex1, float3 vertex2, float3 vvn)
{
  float p0, p1, r;
  float3 m = ccc->m;
  float3 n = ccc->n;
  float3 aabb_center = (m + n)/2;
  float3 extents = (m - n)/2;//m - aabb_center;

  float3 v0 = vertex0 - aabb_center;

  float3 f0 = normalize(vertex1 - vertex0);
  float3 f2 =  - normalize(vertex2 - vertex0);

  float3 a00 = make_float3(0, -f0.z, f0.y);// (1 0 0)xf0
  float3 a02 = make_float3(0, -f2.z, f2.y);// (0 1 0)xf0
  float3 a10 = make_float3(f0.z, 0, -f0.x);// (0 0 1)xf0
  float3 a12 = make_float3(f2.z, 0, -f2.x);// (1 0 0)xf2
  float3 a20 = make_float3(-f0.y, f0.x, 0);// (0 1 0)xf2
  float3 a22 = make_float3(-f2.y, f2.x, 0);// (0 0 1)xf2

  // Test axis a00
  p0 = dot(v0, a00);
  p1 = -dot(f2, a00);
  r = extents.y * fabs(f0.z) + extents.z * fabs(f0.y);

  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }

  // Test axis a02
  p0 = dot(v0, a02);
  p1 =  dot(f0, a02);
  r = extents.y * fabs(f2.z) + extents.z * fabs(f2.y);

  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }

  // Test axis a10
  p0 = dot(v0, a10);
  p1 = -dot(f2, a10);
  r = extents.x * fabs(f0.z) + extents.z * fabs(f0.x);
  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }


  // Test axis a12
  p0 = dot(v0, a12);
  p1 =  dot(f0, a12);
  r = extents.x * fabs(f2.z) + extents.z * fabs(f2.x);

  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }

  // Test axis a20
  p0 = dot(v0, a20);
  p1 = -dot(f2, a20);
  r = extents.x * fabs(f0.y) + extents.y * fabs(f0.x);

  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }


  // Test axis a22
  p0 = dot(v0, a22);
  p1 =  dot(f0, a22);
  r = extents.x * fabs(f2.y) + extents.y * fabs(f2.x);

  if( (p1<0?-p0:p0) > r )
  {
    return false;
  }

  float3 fn = normalize(cross(f2, f0));
  if( fabsf(dot(fn,v0)) > dot(extents,f3abs(fn)) )
    return false;

  //{
  //  float r0 = fmax(extents.x,fmax(extents.y,extents.z))*1.74f * 2;//length(extents); 
  //  float l = length(v0);
  //  if(l>r0)
  //  {
  //    if( dot(v0,(ccc->conec_down))/l > ccc->cos_coner_down )
  //      return false;
  //    if( dot(v0,(ccc->conec_up))/l > ccc->cos_coner_up )
  //      return false;
  //  }
  //}

  float3 v1 = v0 + f0*10;
  float3 v2 = v0 - f2*10;
  if (fmax(v0.x, fmax(v1.x, v2.x)) < -extents.x || fmin(v0.x, fmin(v1.x, v2.x)) > extents.x)
    return false;
  if (fmax(v0.y, fmax(v1.y, v2.y)) < -extents.y || fmin(v0.y, fmin(v1.y, v2.y)) > extents.y)
    return false;
  if (fmax(v0.z, fmax(v1.z, v2.z)) < -extents.z || fmin(v0.z, fmin(v1.z, v2.z)) > extents.z)
    return false;

  return true;
}


__device__ bool ray_aabb_intersect(float3 m, float3 n, curay r)
{
  float3 rpos = r.u0;
  float3 rdir = r.du;
  float t[10];
  t[1] = (n.x - rpos.x) / rdir.x;
  t[2] = (m.x - rpos.x) / rdir.x;
  t[3] = (n.y - rpos.y) / rdir.y;
  t[4] = (m.y - rpos.y) / rdir.y;
  t[5] = (n.z - rpos.z) / rdir.z;
  t[6] = (m.z - rpos.z) / rdir.z;
  t[7] = fmax(fmax(fmin(t[1], t[2]), fmin(t[3], t[4])), fmin(t[5], t[6]));
  t[8] = fmin(fmin(fmax(t[1], t[2]), fmax(t[3], t[4])), fmax(t[5], t[6]));
  t[9] = (t[8] < 0 || t[7] > t[8]) ? 0 : t[7];
  return t[9];
}

__device__ bool ray_triangle_intersect( curay r, float3 v0, float3 v1, float3 v2, float3 &pt )
{
  const float EPSILON = 0.0000001f;
  float3 edge1, edge2, h, s, q;
  float a, f, u, v;
  edge1 = v1 - v0;
  edge2 = v2 - v0;
  h = cross( r.du, edge2);
  a = dot( edge1, h);
  if (a > -EPSILON && a < EPSILON)
    return false;    // This ray is parallel to this triangle.
  f = 1.0 / a;
  s = r.u0 - v0;
  u = f * dot(s,h);
  if (u < 0.0 || u > 1.0)
    return false;
  q = cross(s,edge1);
  v = f * dot(r.du,q);
  if (v < 0.0 || u + v > 1.0)
    return false;
 // At this stage we can compute t to find out where the intersection point is on the line.
  float t = f * dot(edge2,q);
  if (t > EPSILON && t < 1 / EPSILON) // ray intersection
  {
    pt = r.u0 + r.du * t;
   return true;
  }
  else // This means that there is a line intersection but not a ray intersection.
    return false;
}

__device__ bool ray_triangle_intersect( const CUParam &param, int ti, curay r, float3 &pt )
{
  float3 v0 = param.vmap[ param.vimap[3*ti+0] ];
  float3 v1 = param.vmap[ param.vimap[3*ti+1] ];
  float3 v2 = param.vmap[ param.vimap[3*ti+2] ];
  return ray_triangle_intersect( r, v0, v1, v2, pt );
}
//
/////////////////////////////////////////////////////////




__device__ int get_u0_facing( const CUParam &param, const cups &ps, const float3 &pt, int ti, int s, int intersect_ei)
{
  float3 pn0 = normalize(cross(ps.ev0 - ps.origin, ps.ev1 - ps.origin));
  float3 pn1 = normalize(pt - ps.origin);
  float3 pn2 = normalize(cross(pn1, pn0));
  float3 pt1;
  {
    float3 ev0, ev1, ev2;
    ev0 = param.vmap[param.vimap[3 * ti + s]];
    ev1 = param.vmap[param.vimap[3 * ti + (s + 1) % 3]];
    ev2 = param.vmap[param.vimap[3 * ti + (s + 2) % 3]];
    float3 e0 = ev2 - ps.origin;
    float3 e1 = ev1 - ps.origin;
    float3 ee;
    if (dot(e0, pn0) > 0 != dot(e1, pn0) > 0)
      ee = ev1 - ev2;
    else
      ee = ev0 - ev2;
    pt1 = -dot(e0, pn0) / dot(ee, pn0)*ee + ev2;
  }
  bool u0_facing = dot(pt1 - pt, pn2) > 0;
  //return u0_facing;

  int occlusion_count;
  if( param.emap[intersect_ei].t1==-1 )
  {
    if (u0_facing)
      occlusion_count = -1;
    else
      occlusion_count = 1;
  }
  else
  {
    if (u0_facing)
      occlusion_count = -2;
    else
      occlusion_count = 2;
  }

  return occlusion_count;

}

////////////////////////////////////////////////////////
// Line Sample BVH traversal functions
//
__device__ void process_ps_tri_edge_cone(const CUParam &param, cubvh *ccc, const cups &ps
  , EdgeJ *edgej, int &n_edgej
)
{

  int *idx = ccc->idx;
  int ei0 = ps.eidx;

  // for( int i=0; i<ccc->ni; i++ )
  int i=0;
  { 
    int ti = idx[i];
    float3 v0 = param.vmap[ param.vimap[3*ti+0] ];
    float3 n0 = param.nmap[ ti ];

    if (ti != ps.tri_idx)
    {
      int s;
      for (s = 0; s < 3; s++)
      {
        int intersect_ei = param.eimap[3*ti+s];
      
        if( !(dot(n0, ps.origin-v0)>0) && param.emap[intersect_ei].t1!=-1 )
          continue;

        int eivi0, eivi1, ei0vi0, ei0vi1;
        eivi0 = param.emap[intersect_ei].v0;
        eivi1 = param.emap[intersect_ei].v1;
        
        if(ei0!=-1)
        {
          ei0vi0 = param.emap[ei0].v0;
          ei0vi1 = param.emap[ei0].v1;
          if( eivi0 == ei0vi0 || eivi0 == ei0vi1 || eivi1 == ei0vi0 || eivi1 == ei0vi1 )
            continue;
        }

        float3 ev0, en0, en1;
        ev0 = param.vmap[eivi0];
        en0 = param.nmap[param.emap[intersect_ei].t0];

        if(param.emap[intersect_ei].t1!=-1)
          en1 = param.nmap[param.emap[intersect_ei].t1];
        else
          en1 = -en0;

        float3 vv0 = ps.origin - ev0;
        float f0 = dot(vv0, en0);
        float f1 = dot(vv0, en1);

        if( f0>0!=f1>0 )      
        {          
          float es1;
          float3 pt;
          if(
            ps_edge_intersectxa(
              param, ps.vvn, ps.origin, ps.ev0, ps.ev1, ei0, intersect_ei, param.vimap[3*ti+(s+2)%3],
              pt, es1
            )
          ){
            if(n_edgej<INTERSECT_EDGE_SIZE)
            {
              int u0_facing = get_u0_facing(param, ps, pt, ti, s, intersect_ei );
              float es0 = LsPsXRatio(param, intersect_ei, ps.origin, ps.ev0, ps.ev1 );
              es0 = clamp( es0,-.5f,.5f);
              edgej[n_edgej].es0 = es0;
              edgej[n_edgej].u0_facing = u0_facing;
              edgej[n_edgej].es1  = es1;
              edgej[n_edgej].ei1 = intersect_ei;
              n_edgej++;
            }
          }

        }
      }
    }
  }
}

__device__ void ps_bvh_flattened_cone( const CUParam &param, const cups &ps, int addr,
  EdgeJ *edgej, int &n_edgej
){
  cubvh *ccc = (cubvh*) (param.cbvh + addr);
  int my_note[PS_BVH_FLATTENED_CONE_NOTE_SIZE];
  int n = 0;
  my_note[n++] = addr;
  while(n && n<PS_BVH_FLATTENED_CONE_NOTE_SIZE-1 )
  {
    ccc = (cubvh*) (param.cbvh + my_note[n-1]); 
    n--;

    if(ps_aabb_cone_intersect(ccc, ps.origin, ps.ev0, ps.ev1, ps.vvn))
    {
      if(ccc->ni)
        process_ps_tri_edge_cone( param, ccc, ps, edgej, n_edgej );
      else
      {
        if(ccc->left)
          my_note[n++] = ccc->left;
        if(ccc->right)
          my_note[n++] = ccc->right;
      }
    }
  }
}
//
/////////////////////////////////////////////////////////



////////////////////////////////////////////////////////
// Probing Ray BVH traversal functions 
//
__device__ void process_nu0( const CUParam &param, cubvh *ccc, curay r, 
  float3 ps_origin, int ei0, int tri_idx, int &fnu0, int &bnu0 )
{      
  int *idx = ccc->idx;
  for( int i=0; i<ccc->ni; i++ )
  { 
    int ti = idx[i];
     float3 v0 = param.vmap[param.vimap[3 * ti + 0]];
     float3 n0 = param.nmap[ti];
    
    if( ei0==-1 || (param.emap[ei0].t0 != ti && param.emap[ei0].t1 != ti) )
    {
      float3 pp;
      if(ray_triangle_intersect(param, idx[i], r, pp ))
      {
        float d = dot(n0, ps_origin-v0);
        if( d>-flt_epslion )
          fnu0 = fnu0 + 1;
        if( d<+flt_epslion )
          bnu0 = bnu0 + 1;
      }
    }
  }
}

__device__ void rt_cal_nu0_flattened(const CUParam &param, curay r, int addr, 
  float3 ps_origin, int ei0, int tri_idx,  int &fnu0, int &bnu0)
{
  cubvh *ccc = (cubvh*) (param.cbvh + addr);

  int my_note[32];
  int n = 0;

  my_note[n++] = addr;

  while(n && n<31 )
  {
    ccc = (cubvh*) (param.cbvh + my_note[n-1]); 
    n--;

    if(ray_aabb_intersect(ccc->m, ccc->n, r))
    {
      if(ccc->ni)
        process_nu0( param, ccc, r, ps_origin, ei0, tri_idx, fnu0, bnu0);
      if (ccc->left)
        my_note[n++] = ccc->left;
      if (ccc->right)
        my_note[n++] = ccc->right;
    }
  }
}
//
/////////////////////////////////////////////////////////




////////////////////////////////////////////////////////
// PointI functions 
//
//__device__ PointI getPointI( const CUParam &param, int ei, float3 pt )
//{
//  float3 ev0 = param.vmap[ param.emap[ei].v0 ];
//  float3 ev1 = param.vmap[ param.emap[ei].v1 ];
//  PointI pp;
//  pp.ei = ei;
//  pp.es = length(pt-ev0)/length(ev1-ev0)-.5;
//  return pp;
//}
//__device__ PointI getPointI( const CUParam &param, const cups &ps, float theta_hp )
//{
//  PointI pp;
//  pp.ei = ps.eidx;
//  pp.es = theta_hp;
//  if( ps.tri_idx != param.emap[ps.eidx].t0 )
//    pp.es = -pp.es;
//  return pp;
//}
//__device__ PointI getPointI_ev0( const CUParam &param, const cups &ps )
//{
//  return getPointI(param,ps,-.5f);
//}
//__device__ PointI getPointI_ev1( const CUParam &param, const cups &ps )
//{
//  return getPointI(param,ps, .5f);
//}

__device__ float3 eval( const CUParam &param, const PointI &pp, const float3 &origin, const float3 &ax, const float3 &ay, const float3 &az )
{
  if(pp.ei!=-1)
  {
    float3 ev0 = param.vmap[param.emap[pp.ei].v0];
    float3 ev1 = param.vmap[param.emap[pp.ei].v1];
    return pp.es*(ev1-ev0)+(ev1+ev0)/2; 
  }else
  {
    return origin + (cosf(pp.es)*az + sinf(pp.es)*ax)*2;
  }
}


__device__ void extend_memory_record_edge_PointI( int &recount, int n, PointI* &record_edge0, PointI* &record_edge )
{
  if( recount+n>(RECORD_EDGE_SIZE-2) )
  {
    //int pos = atomicAdd( extended, 32 );
    //float3 *tmp = &((float3*)&extended[1])[pos];
    //((int*)record_edge0)[0] = recount;
    //((int*)record_edge0)[1] = pos;

    printf( "record_edge0 memory exceeded\n" );

    //record_edge0 = tmp;
    ((int*)record_edge0)[0] = 0;
    ((int*)record_edge0)[1] = -1;
    record_edge = record_edge0+2;
    recount = 0;
  }
}
__device__ void record_visible_edge(const CUParam &param, const cups &ps, EdgeJ *edgej, int einum, int n, PointI *record_edge00)
{
  PointI *record_edge;
  PointI *record_edge0 = record_edge00;
  ((int*)record_edge0)[0] = 0;
  ((int*)record_edge0)[1] = -1;
  record_edge = record_edge0+2;
  int &recount = *((int*)record_edge00) = 0;

  float3 ax, ay, az;
  PointI u0, u1, u2;
  int i;
  float tu, tv;
  float hp, t0, t1;
  bool A, B, eeflip;

  eeflip = ps.tri_idx != param.emap[ps.eidx].t0;
  A = dot(ps.ev0-ps.origin, ps.vvn)>0;
  B = dot(ps.ev1-ps.origin, ps.vvn)>0;
  hp = LsPlXRatio(ps.origin, ps.vvn, ps.ev0, ps.ev1);
  t0 = tu = (!A && B) ? hp : -.5f;
  t1 = (A && !B) ? hp : .5f;

  for( i=0; i<einum; i++ )
  {
    tv = edgej[i].es0;           
    if(t0<tv)
      break;
    n += edgej[i].u0_facing;
  }

  for(    ; i<einum; i++ )
  {
    tv = edgej[i].es0;           
    if(t1<tv)
      break;
    if( n==0 && fabsf(tv-tu)>.00002f )
    {
      u0.ei = ps.eidx;
      u0.es = eeflip ? -tu : tu;
      u1.ei = ps.eidx;
      u1.es = eeflip ? -tv : tv;
      u2.ei = edgej[i].ei1;
      u2.es = edgej[i].es1;
      extend_memory_record_edge_PointI( recount, 4, record_edge0, record_edge );
      record_edge[recount++] = u0;
      record_edge[recount++] = u1;
      record_edge[recount++] = u1;
      record_edge[recount++] = u2;
    }
    n += edgej[i].u0_facing;
    tu = tv;
  }

  tv = t1;
  if( n==0 && fabsf(tv-tu)>.00002f )
  {
    u0.ei = ps.eidx;
    u0.es = eeflip ? -tu : tu;
    u1.ei = ps.eidx;
    u1.es = eeflip ? -tv : tv;
    extend_memory_record_edge_PointI( recount, 2, record_edge0, record_edge );
    record_edge[recount++] = u0;
    record_edge[recount++] = u1;
    if( A && !B )
    {
      get_coordinate_system(ps, ax, ay, az);
      u2 = getHemiPointI_es( ps, tv, ax,ay,az );
      extend_memory_record_edge_PointI( recount, 2, record_edge0, record_edge );
      record_edge[recount++] = u1;
      record_edge[recount++] = u2;
    }
  }
}
__device__ void record_visible_edge_plane(const CUParam &param, const cups &ps, EdgeJ *edgej, int einum, int n, PointI *record_edge00)
{
  PointI *record_edge;
  PointI *record_edge0 = record_edge00;
  ((int*)record_edge0)[0] = 0;
  ((int*)record_edge0)[1] = -1;
  record_edge = record_edge0+2;
  int &recount = *((int*)record_edge0) = 0;

  float3 ax, ay, az;
  PointI u0, u1, u2;
  int i;
  float tu, tv;

  get_coordinate_system(ps, ax, ay, az);
  tu = -.5f;
  for( i=0; i<einum; i++ )
  {
    tv = edgej[i].es0;
    if( n==0 && fabsf(tv-tu)>.00002f )
    {
      extend_memory_record_edge_PointI( recount, 4, record_edge0, record_edge );
      u0 = getHemiPointI_es( ps,tu, ax,ay,az );
      u1 = getHemiPointI_es( ps,tv, ax,ay,az );
      u2.ei = edgej[i].ei1;
      u2.es = edgej[i].es1;
      record_edge[recount++] = u0;
      record_edge[recount++] = u1;
      record_edge[recount++] = u1;
      record_edge[recount++] = u2;
    }
    n += edgej[i].u0_facing;
    tu = tv;
  }
  tv = .5f;
  if( n==0 && fabsf(tv-tu)>.00002f )
  {
    extend_memory_record_edge_PointI( recount, 2, record_edge0, record_edge );
    u0 = getHemiPointI_es( ps,tu, ax,ay,az );
    u1 = getHemiPointI_es( ps,tv, ax,ay,az );
    record_edge[recount++] = u0;
    record_edge[recount++] = u1;
  }
  
  //*((int*)record_edge0) = recount;
}




__device__ void find_pst_visible_region(const CUParam &param, const cups &ps 
  , PointI *CURecord_vis, EdgeJ *edgej, int &n_edgej)
{
  if(n_edgej>INTERSECT_EDGE_SIZE)
  {
    n_edgej = INTERSECT_EDGE_SIZE;
    printf( "find_pst_visible_region(), buffer overflowed, %s.\n", "theta[]" );
  }
  cg_quick_sort(edgej, n_edgej);

  float thetam;
  {
    float md, d;
    float tu, tv;
    float pu, pv;
    tu = -.5f; //0;
    md = -flt_max;
    for( int i=0; i<n_edgej; i++ )
    {
      tv = edgej[i].es0;
      d  = tv - tu;

      if( md<d )
      {
        md = d;
        pu = tu;
        pv = tv;
      }
      tu = tv;
    }
    tv = .5f;//1;
    d  = tv - tu;
    if( md<d )
    {
      md = d;
      pu = tu;
      pv = tv;
    }
    thetam = (pu+pv)/2.f;
  }
  //thetam = 0;

  int nu0;
  {
    float3 vertexm;
      //vertexm = (ps.ev1-ps.ev0)*(thetam+.5f)+ps.ev0;
      vertexm = es_to_float3( thetam, ps.ev0, ps.ev1 );
      //vertexm = es_to_float3( 0.f, ps.ev0, ps.ev1 );

    curay rnu0;
    int fnu0, bnu0;
      rnu0.u0 = ps.origin;
      rnu0.du = normalize(vertexm - ps.origin);
      fnu0 = 0;
      bnu0 = 0;
    rt_cal_nu0_flattened(param, rnu0, 0, ps.origin, ps.eidx, ps.tri_idx, fnu0, bnu0);
    nu0 = fnu0 + bnu0;
   }

  int n = 0;
  for( int i=0; i<n_edgej && edgej[i].es0<thetam; i++ )
    n += edgej[i].u0_facing;
  n = nu0 - n;

  if( ps.eidx!=-1 )
    record_visible_edge(param, ps, edgej, n_edgej, n, CURecord_vis);
  else
    record_visible_edge_plane(param, ps, edgej, n_edgej, n, CURecord_vis);

}



__global__ void call_ps_bvh_flattened_cone_speed3( CUParam param, cups *ps,  PointI *CURecord_vis, int nps, int ps0)
{
  // gridDim.x 256; [0,(nps+(256-1))/256) <= blockIdx.x;
  // blockDim.x 256; [0,256) <= threadIdx.x;
  int psi = blockIdx.x * blockDim.x + threadIdx.x + ps0;
  if( psi>=nps )
    return;

  EdgeJ edgej[INTERSECT_EDGE_SIZE];
  int n_edgej = 0;

  ps_bvh_flattened_cone(param, ps[psi], 0, edgej, n_edgej);
  find_pst_visible_region(param, ps[psi], &CURecord_vis[psi*RECORD_EDGE_SIZE], edgej, n_edgej);
}


__global__ void tightpack_vis_cal_offset(int *vis_offset, PointI *vis_buf, int nps)
{
  
  // blockDim.x 256;  [0,  32) <- threadIdx.x
  // gridDim.x (nps+(256-1))/256;  [0, (nps+(256-1))/256) <- blockIdx.x

  int id = blockIdx.x*blockDim.x + threadIdx.x;
  if(id>=nps)
    return;

  vis_offset[id] = *((int*)&vis_buf[id*RECORD_EDGE_SIZE]);
}


__global__ void tightpack_vis_copy(PointI *vis_tightpacked, PointI *vis_buf, int *vis_offset, int nps)
{

  // blockDim.x 256;  [0,  32) <- threadIdx.x
  // gridDim.x (nps+(256-1))/256;  [0, (nps+(256-1))/256) <- blockIdx.x
  int id = blockIdx.x*blockDim.x + threadIdx.x;

  if(id>=nps)
    return;

  int ii = vis_offset[id+0];
  int ni = vis_offset[id+1]-vis_offset[id+0];
  PointI *dat = &vis_buf[id*RECORD_EDGE_SIZE]+2;
  for( int i=0; i<ni; i++ )
    vis_tightpacked[ii+i] = dat[i];

}

__global__ void tightpack_vis_cal_vinfo(int *vis_tightpacked_vinfo, int *vis_offset, int *front_back_edge_vinfo, int n_vertex)
{
  // blockDim.x 256;  [0,  32) <- threadIdx.x
  // gridDim.x (n_vertex+1+(256-1))/256;  [0, (n_vertex+1+(256-1))/256) <- blockIdx.x

  int vi = blockIdx.x*blockDim.x + threadIdx.x; 
  if(vi>=n_vertex)
    return;
  vis_tightpacked_vinfo[vi] = vis_offset[   front_back_edge_vinfo[vi] + vi*32   ];
}


__host__ void my_check_cuda()
{
  cudaError_t cudaStatus;
  cudaStatus = cudaGetLastError();
  if( cudaStatus != cudaSuccess )
  {
    printf( "kernel launch error: %s\n", cudaGetErrorString(cudaStatus));
    exit(0);
  }
  cudaStatus = cudaDeviceSynchronize();
  if( cudaStatus != cudaSuccess )
  {
    printf( "cuda sync error: %s\n", cudaGetErrorString(cudaStatus));
    exit(0);
  }
}

__host__ void cu_ps_bvh_flattened_cone_speed_g0(
  PointI *&vis_tightpacked, int &vis_tightpacked_size, int &_vis_tightpacked_size, int *vis_tightpacked_vinfo, 
  PointI *vis_buf, int *vis_offset,
  const CUParam &param, cups *ps_buf, int *front_back_edge_vinfo, int nps,
  int n_vertex
){
  //my_check_cuda();

  // depends on ps0 modify call_ps_bvh_flattened_cone_speed3()
  call_ps_bvh_flattened_cone_speed3<<< (nps+(256-1))/256,256 >>>
    (param, ps_buf, vis_buf, nps, 0);
  //printf( "  call_ps_bvh_flattened_cone_speed3 done\n" );
  //my_check_cuda();


  tightpack_vis_cal_offset<<<  (nps+(256-1))/256, 256  >>>(vis_offset+1, vis_buf, nps);
  //printf( "  tightpack_vis_cal_offset done\n" );
  //my_check_cuda();

  thrust::inclusive_scan(thrust::device, vis_offset+1, vis_offset+1 + nps, vis_offset+1);
  //printf( " CUDA thrust::inclusive_scan done\n" );
  //my_check_cuda();

  cudaMemcpy( &vis_tightpacked_size, vis_offset+nps, sizeof(int), cudaMemcpyDeviceToHost );
  //printf( " CUDA vis_tightpacked_size %i\n", vis_tightpacked_size );
  //my_check_cuda();

  if( _vis_tightpacked_size<vis_tightpacked_size )
  {
    _vis_tightpacked_size = int(vis_tightpacked_size*1.5);
    if(vis_tightpacked)
      cudaFree(vis_tightpacked);

    //printf( "  vis_tightpacked %f MB\n", double(_vis_tightpacked_size)*sizeof(PointI)/1024/1024 );
    cudaMalloc(&vis_tightpacked, _vis_tightpacked_size*sizeof(PointI));
    //my_check_cuda();
  }

  tightpack_vis_copy<<< (nps+(256-1))/256,256 >>>( vis_tightpacked, vis_buf, vis_offset, nps );
  //my_check_cuda();

  tightpack_vis_cal_vinfo<<< (n_vertex+1+(256-1))/256, 256 >>>(vis_tightpacked_vinfo, vis_offset, front_back_edge_vinfo, n_vertex+1);
  //my_check_cuda();
}





__global__ void ps_bvh_flattened_cone_testing1( CUParam param, cups *ps,  PointI *CURecord_vis, int nps, int ps0)
{
  // gridDim.x 256; [0,(nps+(256-1))/256) <= blockIdx.x;
  // blockDim.x 256; [0,256) <= threadIdx.x;
  int psi = blockIdx.x * blockDim.x + threadIdx.x + ps0;
  if( psi>=nps )
    return;

  EdgeJ edgej[INTERSECT_EDGE_SIZE];
  int n_edgej = 0;

  ps_bvh_flattened_cone(param, ps[psi], 0, edgej, n_edgej);
}

__global__ void ps_bvh_flattened_cone_testing2( CUParam param, cups *psbuf,  PointI *CURecord_vis, int nps, int ps0, int *fnu0_buf, int *bnu0_buf)
{
  // gridDim.x 256; [0,(nps+(256-1))/256) <= blockIdx.x;
  // blockDim.x 256; [0,256) <= threadIdx.x;
  int psi = blockIdx.x * blockDim.x + threadIdx.x + ps0;
  if( psi>=nps )
    return;

  cups ps;
  float3 vertexm;
  curay rnu0;
  int fnu0, bnu0;

  ps = psbuf[psi];
  vertexm = (ps.ev0+ps.ev1)/2;
  rnu0.u0 = ps.origin;
  rnu0.du = normalize(vertexm - ps.origin);
  fnu0 = 0;
  bnu0 = 0;

  rt_cal_nu0_flattened(param, rnu0, 0, ps.origin, ps.eidx, ps.tri_idx, fnu0, bnu0);

  fnu0_buf[psi] = fnu0;
  bnu0_buf[psi] = bnu0;
}


#include <time.h>

__host__ void cu_ps_bvh_flattened_cone_speed_testing(
  PointI *vis_buf, int *vis_offset,
  const CUParam &param, cups *ps_buf, int *front_back_edge_vinfo, int nps,
  int n_vertex, int n_update
){

  int vertex0;
  int ps0, nps_update;

  vertex0  = n_vertex-n_update;
  cudaMemcpy( &ps0, front_back_edge_vinfo+vertex0, sizeof(int), cudaMemcpyDeviceToHost );
  ps0 = ps0 + 32*vertex0;
  nps_update = nps - ps0;
  //printf( "%i -> %i\n", nps, nps_update );

  int i;
  int nprocess = 100;

  clock_t t0, t1;
  t0 = clock();
  for( i=0; i<nprocess; i++ )
  {
    // depends on ps0 modify call_ps_bvh_flattened_cone_speed3()
    ps_bvh_flattened_cone_testing1<<< (nps_update+(256-1))/256,256 >>>
      (param, ps_buf, vis_buf, nps, ps0);
    cudaDeviceSynchronize();
  }
  t1 = clock();
  printf( "Line sample BVH traversal: average processing time per pass %f s\n", float(t1-t0)/CLOCKS_PER_SEC/nprocess );

  int *fnu0_buf, *bnu0_buf;
  cudaMalloc( &fnu0_buf, ((n_vertex+(256-1))/256) * 256 * sizeof(int) );
  cudaMalloc( &bnu0_buf, ((n_vertex+(256-1))/256) * 256 * sizeof(int) );

  t0 = clock();
  for( i=0; i<nprocess; i++ )
  {
    ps_bvh_flattened_cone_testing2<<< (nps_update+(256-1))/256,256 >>>
      (param, ps_buf, vis_buf, nps, ps0, fnu0_buf, bnu0_buf);
    cudaDeviceSynchronize();
  }
  t1 = clock();

  printf( "Probing ray BVH traversal: average processing time per pass %f s\n", float(t1-t0)/CLOCKS_PER_SEC/nprocess );
  printf( "Number of line sample %i\n", nps_update);

  cudaFree(fnu0_buf);
  cudaFree(bnu0_buf);
  
}





////////////////////////////////////////////////////////
// Not in used functions
//
__device__ void process_node_triangles( const CUParam &param, cubvh *ccc, curay r, CUType0 *res )
{      
  int *idx = ccc->idx;
  for( int i=0; i<ccc->ni; i++ )
  { 
    float3 pt;
    if( ray_triangle_intersect(param, idx[i], r, pt ) )
    {
      if( 
        res->ti<0 || 
        length2(pt-r.u0)<length2(res->pt-r.u0)
      ){
        res->pt = pt;
        res->ti = idx[i];
      }
    }
  }
}
__device__ void rt_bvh_flattened( const CUParam &param, curay r, int addr, CUType0 *res)
{
  cubvh *ccc = (cubvh*) (param.cbvh + addr);

  int my_note[32];
  int n = 0;

  my_note[n++] = addr;

  while(n && n<31 )
  {
    ccc = (cubvh*) (param.cbvh + my_note[n-1]); 
    n--;

    if(ray_aabb_intersect(ccc->m, ccc->n, r))
    {
      if(ccc->ni)
        process_node_triangles( param, ccc, r, res );
      if (ccc->left)
        my_note[n++] = ccc->left;
      if (ccc->right)
        my_note[n++] = ccc->right;
    }
  }
}
__global__ void call_rt_bvh( CUParam param, curay *rays, CUType0 *results )
{
  int i = threadIdx.x;
  results[i].ti = -1;
  rt_bvh_flattened(param, rays[i], 0, &results[i]);
}
__host__ void cu_rt_bvh( const CUParam &param, curay *rays, CUType0 *results, int ns )
{
  call_rt_bvh<<<1,1>>>( param, rays, results);

  cudaError_t cudaStatus;

  cudaStatus = cudaGetLastError();
  if( cudaStatus != cudaSuccess )
    printf( "kernel launch error: %s\n", cudaGetErrorString(cudaStatus));

  cudaStatus = cudaDeviceSynchronize();
  if( cudaStatus != cudaSuccess )
    printf( "cuda sync error: %s\n", cudaGetErrorString(cudaStatus));
}
////////////////////////////////////////////////////////
