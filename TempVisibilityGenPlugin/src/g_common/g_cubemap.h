#ifndef G_CUBEMAP_H
#define G_CUBEMAP_H

#include "g_vector.h"
#include "g_pfm.h"

const unsigned int CM_GL_FACE[6] = {
  0x8515,  //GL_TEXTURE_CUBE_MAP_POSITIVE_X
  0x8516,  //GL_TEXTURE_CUBE_MAP_NEGATIVE_X
  0x8517,  //GL_TEXTURE_CUBE_MAP_POSITIVE_Y
  0x8518,  //GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
  0x8519,  //GL_TEXTURE_CUBE_MAP_POSITIVE_Z
  0x851A,  //GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
};

enum CUBE_MAP_FACE
{
  POS_X,
  NEG_X,
  POS_Y,
  NEG_Y,
  POS_Z,
  NEG_Z
};

int cm_lookup( FLOAT3 v0, int cm_side );
void cm_lookup( FLOAT3 v0, int cmside, int *idx, float *val );
template <typename T>
T cm_lookup( FLOAT3 v0, const T *dat, int cmside )
{
  int idx[4]; float val[4];
  cm_lookup( v0, cmside, idx, val );
  T col;
    col  = dat[idx[0]]*val[0];
    col += dat[idx[1]]*val[1];
    col += dat[idx[2]]*val[2];
    col += dat[idx[3]]*val[3];
  return col;
}

int sm_lookup( FLOAT3 v0, int smw, int smh );
void sm_lookup( FLOAT3 v0, int smw, int smh, int *idx, float *val );
template <typename T>
T sm_lookup( FLOAT3 v0, const T *dat, int smw, int smh )
{
  int idx[4]; float val[4];
  sm_lookup( v0, smw, smh, idx, val );
  T col;
    col  = dat[idx[0]]*val[0];
    col += dat[idx[1]]*val[1];
    col += dat[idx[2]]*val[2];
    col += dat[idx[3]]*val[3];
  return col;
}

void em_lookup( FLOAT3 v0, int em_side, int *idx, float *val );
template <typename T>
T em_lookup( FLOAT3 v0, const T *dat, int emside )
{
  int idx[4]; float val[4];
  em_lookup( v0, emside, idx, val );
  T col;
    col  = dat[idx[0]]*val[0];
    col += dat[idx[1]]*val[1];
    col += dat[idx[2]]*val[2];
    col += dat[idx[3]]*val[3];
  return col;
}

void cm2dp( const GPfm &cmsrc, GPfm &dmsrc, int dpside );

void cm2sm( const GPfm &cm, GPfm &sm, int smw, int smh );
void cm2sm( const GPf1 &cm, GPf1 &sm, int smw, int smh );

void sm2cm( const GPfm &sm, GPfm &cm, int cmside );
void sm2cm( const GPf1 &sm, GPf1 &cm, int cmside );

void get_cubemap( int cm_side, GPdm &cm_dir );
void get_cubemap( int cm_side, GPfm &cm_dir );
void get_smap( int sm_side, GPfm &sm_dir );

void get_cubemap_solid_angle( int cm_side, GPf1 &sld );

void load_tcm( const char *spath, GPfm &cm_src, int cm_side=0 );
void save_tcm( const char *tcm_path, const GPfm &cm_src, int cm_side=0, const char *format="pfm" );

double get_discrepancy( const FLOAT3 *vec, int n_vec );
float area_spherical_trianglef( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 );
float area_spherical_triangle( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 );
float area_triangle( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 );

int*** edge_cmidx( int cm_side );
void prepare_emsrc( const GPfm &cm_src, GPfm &em_src );

void get_smap( int sm_side, GPdm &sm_dir );

void prepare_emsrc( const GPdm &cmsrc, GPdm &emsrc );

void get_cubemap( int cm_side, FLOAT3 *cm_dir );
void get_cubemap_solid_angle( int cm_side, float *val );
void get_smap( int sm_side, FLOAT3 *sm_dir );
void get_smap_solid_angle( int sm_side, float *val );
void get_smap_solid_angle( int sm_side, GPf1 &sld );

void cm2sm_nearest( const FLOAT3 *cm_data, int cm_side, FLOAT3 *sm_data, int sm_width, int sm_height );
void get_hemi_cubemap( int cm_side, FLOAT3 *hm_dir );
void cm2sm( const FLOAT3 *cm_data, int cm_side, FLOAT3 *sm_data, int sm_width, int sm_height );
void sm2cm( FLOAT3 *cm_data, int cm_side, const FLOAT3 *sm_data, int sm_width, int sm_height );
void cm_save( const char *spath, const FLOAT3 *cm_dat, int cm_side );

void get_hammersley( int n_point, GPfm &hsdir );

#endif
