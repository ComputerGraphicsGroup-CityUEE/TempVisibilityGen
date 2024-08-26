#ifndef MY_RAYTRACE_DATATYPES_H
#define MY_RAYTRACE_DATATYPES_H

#include "glm/glm.hpp"


//struct RayPayload
//{
//  vec3 color;
//  float distance;
//
//  vec3 normal; 
//  float reflector;
//
//  vec3 bc;
//  int ti;
//
//  vec3 fn;
//  uint ray_type;
//  vec3 kd;
//  float roughness;
//  vec3 ks;
//  float shininess;
//
//  uint hitkind;
//  uint mtli;
//  vec2 dummy0;
//};
//
//
//struct Material
//{
//  vec4 baseColorFactor;
//  int baseColorTexture;
//  float metallic;
//  float roughness;
//
//  float kd;
//  float ks;
//  float shininess;
//  float cone_radius;
//
//  float dummy02;
//};
//
//struct AnnenLit 
//{
//  mat4 viewInverse;
//  mat4 projInverse;
//  vec3 le;
//  float lightsize;
//  
//  mat4 viewM;
//  mat4 projM;
//  vec3 lc;
//  float dummy0;
//  //mat3 viewM;
//  //float dummy0;
//  //mat3 projM;
//  //float dummy1;
//  //vec3 lc;
//  //float dummy2;
//};
//
//struct CameraProperties 
//{
//  mat4 viewInverse;
//  mat4 projInverse;
//  mat4 mvp;
//  mat4 mvp_prev;
//
//  vec3 epos;
//  int ray_per_pixel;
//
//  uint frame_id;
//  uint mode;
//  int n_bounce;
//  float atteun_c1;
//
//  float atteun_c2;
//  uint bitterli_mode;
//  int bitterli_k;
//  float mis_clight;
//
//  float annen_fs0;
//  float annen_fs1;
//  float annen_fs2;
//  float annen_fs3;
//
//  float render_tone_val;
//  uint mtl_bitmode;
//  float dummy01;
//  float dummy02;
//
//};

struct VisInfo
{
  int eidx;
  float dis;
};

//struct GEdge
//{
//  int v0, t0;
//  int v1, t1;
//};

struct GEdge
{
  int v0, v1;
  int t0, t1;
};


struct EdgeI
{
  int   ei0;
  float es0;
  int   ei1;
  float es1;
};

struct MyVertex
{
  //glm::vec3 pos;
  //glm::vec3 normal;

  glm::vec3 pos;
  float dummy0;
  glm::vec3 normal;
  float dummy1;
  glm::vec2 uv0;
  glm::vec2 uv1;
  glm::vec3 color0;
  float dummy2;


  //glm::vec3 v0;
  //glm::vec3 n0;

  //float dummy0;
  //float dummy1;
  //vec2 uv0;
  //vec2 uv1;
  //vec3 color0;
  //float dummy2;
};

//struct MyPushconst 
//{
//  uint frame_seed;
//  int accid;
//};
//

#endif
