#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>


#include "g_common.h"


#include "g_cubemap.h"


void cm2sm( const GPfm &cm, GPfm &sm, int smw, int smh )
{
  float theta, phi;
  FLOAT3 l;
  int i,j;
  sm.load( smw, smh );
  for( j=0; j<sm.h; j++ )
    for( i=0; i<sm.w; i++ )
    {
      theta = (j+.5f) *     G_PI / sm.h;
      phi   = (i+.5f) * 2 * G_PI / sm.w;
      Angle2Float3( theta, phi, &l );
      sm.pm[j][i] = cm_lookup( l, cm.fm, cm.w );
    }
}

void cm2sm( const GPf1 &cm, GPf1 &sm, int smw, int smh )
{
  float theta, phi;
  FLOAT3 l;
  int i,j;
  sm.load( smw, smh );
  for( j=0; j<sm.h; j++ )
    for( i=0; i<sm.w; i++ )
    {
      theta = (j+.5f) *     G_PI / sm.h;
      phi   = (i+.5f) * 2 * G_PI / sm.w;
      Angle2Float3( theta, phi, &l );
      sm.pm[j][i] = cm_lookup( l, cm.fm, cm.w );
    }
}

void sm2cm( const GPfm &sm, GPfm &cm, int cmside )
{
  GPfm cmdir;
  int i, l;
  cm.load( cmside, 6*cmside );
  get_cubemap( cmside, cmdir );
  l = cmside * cmside * 6;
  for( i=0; i<l; i++ )
    cm.fm[i] = sm_lookup( cmdir.fm[i], sm.fm, sm.w, sm.h );
}

void sm2cm( const GPf1 &sm, GPf1 &cm, int cmside )
{
  GPfm cmdir;
  int i, l;
  cm.load( cmside, 6*cmside );
  get_cubemap( cmside, cmdir );
  l = cmside * cmside * 6;
  for( i=0; i<l; i++ )
    cm.fm[i] = sm_lookup( cmdir.fm[i], sm.fm, sm.w, sm.h );
}

void get_cubemap( int cm_side, GPfm &cm_dir )
{
  cm_dir.load( cm_side, 6*cm_side );
  get_cubemap( cm_side, cm_dir.fm );
}

void get_smap( int sm_side, GPfm &sm_dir )
{
  sm_dir.load( sm_side*2, sm_side );
  get_smap( sm_side, sm_dir.fm );
}

void get_smap( int sm_side, GPdm &sm_dir )
{
  int i,j;
  double theta, phi;

  sm_dir.load( sm_side*2, sm_side );
  for( j=0; j<sm_dir.h; j++ )
  for( i=0; i<sm_dir.w; i++ )
  {
    theta = (j+.5) *     D_PI / sm_dir.h;
    phi   = (i+.5) * 2 * D_PI / sm_dir.w;
    Angle2Double3( theta, phi, &sm_dir.pm[j][i] );
  }
}



void get_cubemap_solid_angle( int cm_side, GPf1 &sld )
{
  sld.load( cm_side, cm_side*6 );
  get_cubemap_solid_angle( cm_side, sld.fm );
}

void load_tcm( const char *spath, GPfm &cm_src, int cm_side )
{
  char face_tag[6][6] = 
  {
    "POS_X", "NEG_X",
    "POS_Y", "NEG_Y",
    "POS_Z", "NEG_Z",
  };

  GPath gp = parse_spath( spath );

  char fpath[6][256];
  char str[256], *tstr;
  int i;
  int cmtype = 0;

  FILE *f0 = fopen( spath, "rt" );
    if( f0==0 )
    {
      printf( "[Error] load_tcm(), file %s not found.\n" ,spath );
      exit(-1);
    }
    while( fgets( str, 256, f0 ) )
    {
      replace_char( str, '\r', 0 );
      replace_char( str, '\n', 0 );
      tstr = strtok( str, " " );
      if (strcmp(tstr, "CM_TYPE") == 0)
      {
        if (strcmp(strtok(NULL, " "), "NVCM") == 0)
          cmtype = 1;
      }else
      {
        for (i = 0; i < 6; i++)
          if (strcmp(tstr, face_tag[i]) == 0)
          {
            sprintf(fpath[i], "%s%s", gp.dname, strtok(NULL, " "));
            break;
          }
      }
    }
  fclose(f0);

  GPfm blk;
    blk.load( fpath[0] );
    if( cm_side==0 )
      cm_side = blk.w;

  cm_src.load( cm_side, cm_side*6 );
  for( i=0; i<6; i++ )
  {
    blk.load( fpath[i] );
    blk.resample( cm_side, cm_side );
    cm_src.draw( blk, 0, i*cm_side );
  }
  if (cmtype == 1)
    cm_src.flip_horizontal();
}

void save_tcm( const char *tcm_path, const GPfm &cm_src, int cm_side, const char *format )
{
  GPath gp = parse_spath( tcm_path );
  char sface[6][3] = { "xp", "xn", "yp", "yn", "zp", "zn" };
  char tface[6][6] = { "POS_X", "NEG_X", "POS_Y", "NEG_Y", "POS_Z", "NEG_Z" };
  char spath[256];
  GPfm blk;
  int i;

  for( i=0; i<6; i++ )
  {
    sprintf( spath, "%s%s.%s.%s", gp.dname, gp.fname, sface[i], format );
    cm_src.getblk( blk, 0, i*cm_src.w, cm_src.w, cm_src.w );
    if( cm_side )
      blk.resample( cm_side, cm_side );
    blk.save( spath, format );
  }

  FILE *f0 = fopen( tcm_path, "wt" );
    fprintf( f0, "CM_TYPE TTCM\n" );
    for( i=0; i<6; i++ )
      fprintf( f0, "%s %s.%s.%s\n", tface[i], gp.fname, sface[i], format );
  fclose(f0);
}

double get_discrepancy( const FLOAT3 *vec, int n_vec )
{
  double beta;
  int i, j;

  beta=0;
  for( j=0; j<n_vec; j++ )
  for( i=0; i<n_vec; i++ )
    beta += 1-2*log(1+sqrt((1.0000015-vdot(vec[j],vec[i]))*0.5));
  beta = sqrt(beta)/(2*sqrt(G_PI)*n_vec);

  return beta;
}

float area_spherical_trianglef( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 )
{
  float t0, t1, t2, area;

  FLOAT3 n0, n1, n2;

  n0 = vcross( p0,p1 );
  n1 = vcross( p1,p2 );
  n2 = vcross( p2,p0 );
  vnormalize( &n0 );
  vnormalize( &n1 );
  vnormalize( &n2 );

  t0 = acosf( vdot( -n0, n2 ) );
  t1 = acosf( vdot( -n1, n0 ) );
  t2 = acosf( vdot( -n2, n1 ) );

  area = t0 + t1 + t2 - G_PI;

  return area;
}

float area_spherical_triangle( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 )
{
  double t0, t1, t2, area;
  double n0[3], n1[3], n2[3];
  double d0, d1, d2;
  FLOAT3 a, b;

  a=p0; b=p1;
  n0[0] =  a.y*b.z - b.y*a.z;
  n0[1] = -a.x*b.z + b.x*a.z;
  n0[2] =  a.x*b.y - b.x*a.y;
  a=p1; b=p2;
  n1[0] =  a.y*b.z - b.y*a.z;
  n1[1] = -a.x*b.z + b.x*a.z;
  n1[2] =  a.x*b.y - b.x*a.y;
  a=p2; b=p0;
  n2[0] =  a.y*b.z - b.y*a.z;
  n2[1] = -a.x*b.z + b.x*a.z;
  n2[2] =  a.x*b.y - b.x*a.y;

  d0=sqrt(n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2]);
  d1=sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
  d2=sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
  if(d0!=0){ n0[0]/=d0; n0[1]/=d0; n0[2]/=d0; }
  if(d1!=0){ n1[0]/=d1; n1[1]/=d1; n1[2]/=d1; }
  if(d2!=0){ n2[0]/=d2; n2[1]/=d2; n2[2]/=d2; }

  d0 = -n0[0]*n2[0]-n0[1]*n2[1]-n0[2]*n2[2];
  d1 = -n1[0]*n0[0]-n1[1]*n0[1]-n1[2]*n0[2];
  d2 = -n2[0]*n1[0]-n2[1]*n1[1]-n2[2]*n1[2];

  d0 = g_clamp( d0, -1.0, 1.0 );
  d1 = g_clamp( d1, -1.0, 1.0 );
  d2 = g_clamp( d2, -1.0, 1.0 );

  t0 = acos( d0 );
  t1 = acos( d1 );
  t2 = acos( d2 );

  area = t0 + t1 + t2 - G_PI;

  return (float)area;
}

float area_triangle( const FLOAT3 &p0, const FLOAT3 &p1, const FLOAT3 &p2 )
{
  float area;
  FLOAT3 l0, l1;
  FLOAT3 a0;

  l0 = p1 - p0;
  l1 = p2 - p0;
  a0 = vcross( l0, l1 );

  area = sqrtf(vdot(&a0,&a0))/2;

  return area;
}

int*** edge_cmidx( int cm_side )
{
  int em_side = cm_side + 2;
  int*** em2cm_idx = (int***) malloc3d( em_side, em_side, 6, sizeof(int) );

  int i, j, k, e;

  for( k=0; k<6; k++ )
    for( j=0; j<em_side; j++ )
      for( i=0; i<em_side; i++ )
        em2cm_idx[k][j][i] = -1;


  for( k=0; k<6; k++ )
    for( j=0; j<cm_side; j++ )
      for( i=0; i<cm_side; i++ )
        em2cm_idx[k][j+1][i+1] = k*cm_side*cm_side + j*cm_side + i;

  int rot[6][4] = {
    {  90, -90,   0,   0 },
    { -90,  90,   0,   0 },
    { 180,   0, -90,  90 },
    {   0, 180,  90, -90 },
    {   0,   0,   0,   0 },
    { 180, 180,   0,   0 },
  };

  int n = em_side-1;
  int efxy[6][4][3] = {
    { {2,0,1}, {3,0,1}, {5,n,1}, {4,0,1} },
    { {2,n,1}, {3,n,1}, {4,n,1}, {5,0,1} },
    { {5,1,0}, {4,1,0}, {0,1,0}, {1,1,0} },
    { {4,1,n}, {5,1,n}, {0,1,n}, {1,1,n} },
    { {2,1,n}, {3,1,0}, {0,n,1}, {1,0,1} },
    { {2,1,0}, {3,1,n}, {1,n,1}, {0,0,1} },
  };

  int m = cm_side-1;
  int exy[4][2][2] = 
  {
    { {0,0}, {m,0} },
    { {0,m}, {m,m} },
    { {0,0}, {0,m} },
    { {m,0}, {m,m} },
  };

  for( k=0; k<6; k++ )
  {
    for( e=0; e<4; e++ )
    {
      switch( rot[k][e] )
      {
        case 0:
          for( j=exy[e][0][1]; j<=exy[e][1][1]; j++ )
            for( i=exy[e][0][0]; i<=exy[e][1][0]; i++ )
              em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+j-exy[e][0][1] ][ efxy[k][e][1]+i-exy[e][0][0] ] = k*cm_side*cm_side + j*cm_side + i;
          break;
        case 180:
          for( j=exy[e][0][1]; j<=exy[e][1][1]; j++ )
            for( i=exy[e][0][0]; i<=exy[e][1][0]; i++ )
              em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+j-exy[e][0][1] ][ efxy[k][e][1]+i-exy[e][0][0] ] = k*cm_side*cm_side + j*cm_side + cm_side-i-1;
          break;
        case 90:
          for( j=exy[e][0][1]; j<=exy[e][1][1]; j++ )
            for( i=exy[e][0][0]; i<=exy[e][1][0]; i++ )
              if( e<2 )
                em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+i-exy[e][0][0] ][ efxy[k][e][1]+j-exy[e][0][1] ] = k*cm_side*cm_side + j*cm_side + i;
              else
                em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+i-exy[e][0][0] ][ efxy[k][e][1]+j-exy[e][0][1] ] = k*cm_side*cm_side + (cm_side-j-1)*cm_side + i;
          break;
        case -90:
          for( j=exy[e][0][1]; j<=exy[e][1][1]; j++ )
            for( i=exy[e][0][0]; i<=exy[e][1][0]; i++ )
            {
              if( e>=2 )
                em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+i-exy[e][0][0] ][ efxy[k][e][1]+j-exy[e][0][1] ] = k*cm_side*cm_side + j*cm_side + i;
              else
                em2cm_idx[ efxy[k][e][0] ][ efxy[k][e][2]+i-exy[e][0][0] ][ efxy[k][e][1]+j-exy[e][0][1] ] = k*cm_side*cm_side + j*cm_side + cm_side-i-1;
            }
          break;
      };
    }
  }

  for( k=0; k<6; k++ )
  {
    em2cm_idx[k][0][0] = em2cm_idx[k][0][1];
    em2cm_idx[k][0][n] = em2cm_idx[k][1][n];
    em2cm_idx[k][n][n] = em2cm_idx[k][n][n-1];
    em2cm_idx[k][n][0] = em2cm_idx[k][n-1][0];
  }

  //for( k=0; k<6; k++ )
  //{
  //  for( j=0; j<em_side; j++ )
  //  {
  //    for( i=0; i<em_side; i++ )
  //      printf( "%3i ", em2cm_idx[k][j][i] );
  //    printf( "\n" );
  //  }
  //  printf( "\n" );
  //}

  return em2cm_idx;
}

void prepare_emsrc( const GPfm &cmsrc, GPfm &emsrc )
{
  int cmside, emside;
  
  cmside = cmsrc.w;
  emside = cmside+2;
  emsrc.load( emside, emside*6 );
  {
    int ***em2cm_idx;
    int i, j, k;
    em2cm_idx = edge_cmidx( cmside );
    for( k=0; k<6; k++ )
    {
      for( j=0; j<emside; j++ )
        for( i=0; i<emside; i++ )
          emsrc.pm[k*emside+j][i] = cmsrc.fm[ em2cm_idx[k][j][i] ];
    }
        //emsrc.flip_horizontal();
        
    FLOAT3 col[8] = {0};
    int vidx[6][4] = 
    {
      { 0, 3, 7, 4 },  // +x face
      { 2, 1, 5, 6 },  // -x face
      { 0, 1, 2, 3 },  // +y face
      { 7, 6, 5, 4 },  // -y face
      { 3, 2, 6, 7 },  // +z face
      { 1, 0, 4, 5 },  // -z face
    };

    int pidx[4][2] = { {0,0}, {1,0}, {1,1}, {0,1} };

    for( j=0; j<6; j++ )
    for( i=0; i<4; i++ )
    {
      col[ vidx[j][i] ] += cmsrc.pm[ j*cmside + pidx[i][1]*(cmside-1) ][ pidx[i][0]*(cmside-1) ]/3;
    }

    for( j=0; j<6; j++ )
    for( i=0; i<4; i++ )
    {
      emsrc.pm[ j*emside + pidx[i][1]*(emside-1) ][ pidx[i][0]*(emside-1) ] = col[ vidx[j][i] ];
    }

    free(em2cm_idx);
  }
}

void prepare_emsrc( const GPdm &cmsrc, GPdm &emsrc )
{
  int cmside, emside;
  
  cmside = cmsrc.w;
  emside = cmside+2;
  emsrc.load( emside, emside*6 );
  {
    int ***em2cm_idx;
    int i, j, k;
    em2cm_idx = edge_cmidx( cmside );
    for( k=0; k<6; k++ )
    {
      for( j=0; j<emside; j++ )
        for( i=0; i<emside; i++ )
          emsrc.pm[k*emside+j][i] = cmsrc.fm[ em2cm_idx[k][j][i] ];
    }
        //emsrc.flip_horizontal();
        
    DOUBLE3 col[8] = {0};
    int vidx[6][4] = 
    {
      { 0, 3, 7, 4 },  // +x face
      { 2, 1, 5, 6 },  // -x face
      { 0, 1, 2, 3 },  // +y face
      { 7, 6, 5, 4 },  // -y face
      { 3, 2, 6, 7 },  // +z face
      { 1, 0, 4, 5 },  // -z face
    };

    int pidx[4][2] = { {0,0}, {1,0}, {1,1}, {0,1} };

    for( j=0; j<6; j++ )
    for( i=0; i<4; i++ )
    {
      col[ vidx[j][i] ] += cmsrc.pm[ j*cmside + pidx[i][1]*(cmside-1) ][ pidx[i][0]*(cmside-1) ]/3;
    }

    for( j=0; j<6; j++ )
    for( i=0; i<4; i++ )
    {
      emsrc.pm[ j*emside + pidx[i][1]*(emside-1) ][ pidx[i][0]*(emside-1) ] = col[ vidx[j][i] ];
    }

    free(em2cm_idx);
  }
}




void get_cubemap( int cm_side, GPdm &cm_dir )
{
  DOUBLE3 v[8] = 
  {
    DOUBLE3(  1,  1, -1 ),  // v0
    DOUBLE3( -1,  1, -1 ),  // v1
    DOUBLE3( -1,  1,  1 ),  // v2
    DOUBLE3(  1,  1,  1 ),  // v3
    DOUBLE3(  1, -1, -1 ),  // v4
    DOUBLE3( -1, -1, -1 ),  // v5
    DOUBLE3( -1, -1,  1 ),  // v6
    DOUBLE3(  1, -1,  1 ),  // v7
  };

  int vidx[6][4] = 
  {
    { 0, 3, 7, 4 },  // +x face
    { 2, 1, 5, 6 },  // -x face
    { 0, 1, 2, 3 },  // +y face
    { 7, 6, 5, 4 },  // -y face
    { 3, 2, 6, 7 },  // +z face
    { 1, 0, 4, 5 },  // -z face
  };
  

  DOUBLE3 o,x,y;
  DOUBLE3 dx,dy;
  int i,j,k;

  cm_dir.load( cm_side, 6*cm_side );
  DOUBLE3 *tcm = cm_dir.fm;

  for( k=0; k<6; k++ )
  {
    o = v[vidx[k][0]];
    x = v[vidx[k][1]];
    y = v[vidx[k][3]];

    dx = (x-o)/double(cm_side);
    dy = (y-o)/double(cm_side);
    for( j=0; j<cm_side; j++ )
      for( i=0; i<cm_side; i++, tcm++ )
      {
        *tcm = vnormalize( o + dx*(i+.5f) + dy*(j+.5f)  );
      }
  }
}



void get_cubemap( int cm_side, FLOAT3 *cm_dir )
{
  FLOAT3 v[8] = 
  {
    FLOAT3(  1,  1, -1 ),  // v0
    FLOAT3( -1,  1, -1 ),  // v1
    FLOAT3( -1,  1,  1 ),  // v2
    FLOAT3(  1,  1,  1 ),  // v3
    FLOAT3(  1, -1, -1 ),  // v4
    FLOAT3( -1, -1, -1 ),  // v5
    FLOAT3( -1, -1,  1 ),  // v6
    FLOAT3(  1, -1,  1 ),  // v7
  };

  int vidx[6][4] = 
  {
    { 0, 3, 7, 4 },  // +x face
    { 2, 1, 5, 6 },  // -x face
    { 0, 1, 2, 3 },  // +y face
    { 7, 6, 5, 4 },  // -y face
    { 3, 2, 6, 7 },  // +z face
    { 1, 0, 4, 5 },  // -z face
  };
  

  FLOAT3 o,x,y;
  FLOAT3 dx,dy;
  FLOAT3 Sij;
  int i,j,k;
  FLOAT3 *tcm = cm_dir;

  for( k=0; k<6; k++ )
  {
    o = v[vidx[k][0]];
    x = v[vidx[k][1]];
    y = v[vidx[k][3]];

    dx = (x-o)/float(cm_side);
    dy = (y-o)/float(cm_side);
    for( j=0; j<cm_side; j++ )
      for( i=0; i<cm_side; i++, tcm++ )
      {
        Sij = o + dx*(i+.5f) + dy*(j+.5f);
        vnormalize( &Sij );
        *tcm = Sij;
      }
  }
}

void get_cubemap_solid_angle( int cm_side, float *val )
{
  FLOAT3 v[8] = 
  {
    FLOAT3(  1,  1, -1 ),  // v0
    FLOAT3( -1,  1, -1 ),  // v1
    FLOAT3( -1,  1,  1 ),  // v2
    FLOAT3(  1,  1,  1 ),  // v3
    FLOAT3(  1, -1, -1 ),  // v4
    FLOAT3( -1, -1, -1 ),  // v5
    FLOAT3( -1, -1,  1 ),  // v6
    FLOAT3(  1, -1,  1 ),  // v7
  };

  int vidx[6][4] = 
  {
    { 0, 3, 7, 4 },  // +x face
    { 2, 1, 5, 6 },  // -x face
    { 0, 1, 2, 3 },  // +y face
    { 7, 6, 5, 4 },  // -y face
    { 3, 2, 6, 7 },  // +z face
    { 1, 0, 4, 5 },  // -z face
  };
  

  FLOAT3 o,x,y;
  FLOAT3 dx,dy;
  FLOAT3 p0,p1,p2,p3;
  float w0,w1,w;

  int i,j,k;
  int half_n;
  half_n = (cm_side+1)/2;

  //for( k=0; k<6; k++ )
  //{
    o = v[vidx[0][0]];
    x = v[vidx[0][1]];
    y = v[vidx[0][3]];

    dx = (x-o)/float(cm_side);
    dy = (y-o)/float(cm_side);
    for( j=0; j<half_n; j++ )
      for( i=0; i<half_n; i++ )
      {
        p0 = o + dx*float(i) + dy*float(j);
        p1 = p0 + dx;
        p2 = p0 + dy;
        p3 = p0 + dx + dy;

        vnormalize( &p0 );
        vnormalize( &p1 );
        vnormalize( &p2 );
        vnormalize( &p3 );

        if( cm_side < 512 ) //53 )
        {
          w0 = area_spherical_triangle( p0, p1, p2 );
          w1 = area_spherical_triangle( p3, p1, p2 );
          w = w0 + w1;
        }else
        {
          w0 = area_triangle( p0, p1, p2 );
          w1 = area_triangle( p3, p1, p2 );
          w = w0 + w1;
        }

        val[         j    * cm_side +         i    ] = w;
        val[         j    * cm_side + (cm_side-i-1) ] = w;
        val[ (cm_side-j-1) * cm_side +         i    ] = w;
        val[ (cm_side-j-1) * cm_side + (cm_side-i-1) ] = w;
      }
  //}

  for( k=1; k<6; k++ )
    memcpy( &val[k*cm_side*cm_side], val, cm_side*cm_side*sizeof(float) );
}

void get_smap( int sm_side, FLOAT3 *sm_dir )
{
  int sm_width, sm_height;
    sm_width  = sm_side * 2;
    sm_height = sm_side;

  int i,j;

  float theta, phi;
  FLOAT3 l;
  for( j=0; j<sm_height; j++ )
    for( i=0; i<sm_width; i++ )
    {
      theta = (j+.5f) *     G_PI / sm_height;
      phi   = (i+.5f) * 2 * G_PI / sm_width;
      Angle2Float3( theta, phi, sm_dir++ );
    }
}

void get_smap_solid_angle( int sm_side, float *val )
{
  int sm_width, sm_height;
    sm_width  = sm_side * 2;
    sm_height = sm_side;

  int i,j;
  float sint;
  for( j=0; j<sm_height; j++ )
  {
    sint = sinf( (j+.5f) *     G_PI / sm_height );
    for( i=0; i<sm_width; i++, val++ )
      *val = 2*G_PI*sint/sm_width * G_PI/sm_height;
  }
}

void get_smap_solid_angle( int sm_side, GPf1 &sld )
{
  sld.load( sm_side*2, sm_side );
  get_smap_solid_angle( sm_side, sld.fm );
}






void get_hemi_cubemap( int cm_side, FLOAT3 *hm_dir )
{
  FLOAT3 *cm_dir = (FLOAT3*) malloc( cm_side * cm_side * 6 * sizeof(FLOAT3) );
  get_cubemap( cm_side, cm_dir );

  memcpy( &hm_dir[ 0*cm_side*cm_side/2 ], &cm_dir[ 0*cm_side*cm_side ], cm_side*cm_side/2*sizeof(FLOAT3) );
  memcpy( &hm_dir[ 1*cm_side*cm_side/2 ], &cm_dir[ 1*cm_side*cm_side ], cm_side*cm_side/2*sizeof(FLOAT3) );
  memcpy( &hm_dir[ 2*cm_side*cm_side/2 ], &cm_dir[ 2*cm_side*cm_side ], cm_side*cm_side  *sizeof(FLOAT3) );
  memcpy( &hm_dir[ 4*cm_side*cm_side/2 ], &cm_dir[ 4*cm_side*cm_side ], cm_side*cm_side/2*sizeof(FLOAT3) );
  memcpy( &hm_dir[ 5*cm_side*cm_side/2 ], &cm_dir[ 5*cm_side*cm_side ], cm_side*cm_side/2*sizeof(FLOAT3) );

  free( cm_dir );
}




void cm_save( const char *spath, const FLOAT3 *cm_dat, int cm_side )
{
  int n = cm_side;
  int i;

  const FLOAT3 *face[6];
    for( i=0; i<6; i++ )
      face[i] = cm_dat + i*n*n;

  GPfm tmp;
    tmp.load( 4*n, 3*n );

    tmp.draw( face[0],   0,   n, n,n );
    tmp.draw( face[1], 2*n,   n, n,n );
    tmp.draw( face[2],   n,   0, n,n );
    tmp.draw( face[3],   n, 2*n, n,n );
    tmp.draw( face[4],   n,   n, n,n );
    tmp.draw( face[5], 3*n,   n, n,n );

  tmp.save( spath );
}


void cm2sm( const FLOAT3 *cm_data, int cm_side, FLOAT3 *sm_data, int sm_width, int sm_height )
{
  float theta, phi;
  FLOAT3 l;
  int i,j;
  for( j=0; j<sm_height; j++ )
    for( i=0; i<sm_width; i++ )
    {
      theta = (j+.5f) *     G_PI / sm_height;
      phi   = (i+.5f) * 2 * G_PI / sm_width;
      Angle2Float3( theta, phi, &l );

      sm_data[j*sm_width+i] = cm_lookup( l, cm_data, cm_side );
      //sm_data[j*sm_width+i] = cm_data[cm_lookup( l, cm_side )];
    }
}

void cm2sm_nearest( const FLOAT3 *cm_data, int cm_side, FLOAT3 *sm_data, int sm_width, int sm_height )
{
  float theta, phi;
  FLOAT3 l;
  int i,j;
  for( j=0; j<sm_height; j++ )
    for( i=0; i<sm_width; i++ )
    {
      theta = (j+.5f) *     G_PI / sm_height;
      phi   = (i+.5f) * 2 * G_PI / sm_width;
      Angle2Float3( theta, phi, &l );

      int idx = cm_lookup( l, cm_side );
      sm_data[j*sm_width+i] = cm_data[idx];
    }
}


// FLOAT3 sm_lookup_nearest( FLOAT3 v0, const FLOAT3 *sm_data, int sm_width, int sm_height )
// {
//   float theta, phi;
//   Float32Angle( &theta, &phi, &v0 );
//   if( phi<0 ) phi+= 2*G_PI;
// 
//   float sx, sy;
//   sx = phi / ( 2 * G_PI );
//   sy = theta / G_PI;
// 
//   sx = G_CLAMP( sx, 0, 1-FLT_EPSILON ) * sm_width  + .5f;
//   sy = G_CLAMP( sy, 0, 1-FLT_EPSILON ) * sm_height + .5f;
// 
//   int x0,y0;
//     x0 = (int)G_CLAMP( floor(sx), 1, sm_width  ) - 1;
//     y0 = (int)G_CLAMP( floor(sy), 1, sm_height ) - 1;
// 
//   int idx0;
//     idx0 = y0*sm_width + x0;
// 
//   FLOAT3 val;
//     val = sm_data[idx0];
// 
//   return val;
// }



void sm2cm( FLOAT3 *cm_data, int cm_side, const FLOAT3 *sm_data, int sm_width, int sm_height )
{
  GPfm cm_dir;
    cm_dir.load( cm_side, 6*cm_side );
    get_cubemap( cm_side, cm_dir.fm );

  int i,l;
  l = cm_side * cm_side * 6;

  for( i=0; i<l; i++ )
    cm_data[i] = sm_lookup( cm_dir.fm[i], sm_data, sm_width, sm_height );
}










int cm_lookup( FLOAT3 v0, int cm_side )
{
  FLOAT3 v = v0;

  float n_side = cm_side-FLT_EPSILON;

  int idx;
  int plane,x,y;
  float sx, sy;

  int major_plane =0;
  if( fabs(v.x) < fabs(v.y) )
    major_plane = 1;

  if( major_plane )
  {
    if( fabs(v.y) < fabs(v.z) )
      major_plane = 2;
  }else
  {
    if( fabs(v.x) < fabs(v.z) )
      major_plane = 2;
  }


  switch( major_plane )
  {
    case 0:
      if( v.x > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.x;
        sx = (1 + v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 0;
      }else
      {
        //{ -1,  1,  1 },  // v2
        v = v/(-v.x);
        sx = (1 - v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 1;
      }
      break;

    case 1:
      if( v.y > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.y;
        sx = (1 - v.x) / 2;
        sy = (1 + v.z) / 2;
        plane = 2;
      }else
      {
        //{  1, -1,  1 },  // v7
        v = v/(-v.y);
        sx = (1 - v.x) / 2;
        sy = (1 - v.z) / 2;
        plane = 3;
      }
      break;

    case 2:
      if( v.z > 0 )
      {
        //{  1,  1,  1 },  // v3
        v = v/v.z;
        sx = (1 - v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 4;
      }else
      {
        //{ -1,  1, -1 },  // v1
        v = v/(-v.z);
        sx = (1 + v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 5;
      }
      break;
  };


  sx = G_CLAMP( sx, 0, 1-FLT_EPSILON );
  sy = G_CLAMP( sy, 0, 1-FLT_EPSILON );

  x = (int) (sx * cm_side);
  y = (int) (sy * cm_side);


  idx = cm_side*cm_side*plane + y*cm_side + x;

  if( idx>=cm_side*cm_side*6 )
  {
    printf( "index out of bound... " );
    FILE *f0 = fopen( "log.txt", "at" );
      fprintf( f0, "[Error] cm_lookup(), index out of bound\n" );
      fprintf( f0, "  (x, y, z)     = (%f, %f, %f)\n", v0.x, v0.y, v0.z );
      fprintf( f0, "  (plane, x, y) = (%i, %i, %i)\n", plane, x, y );
    fclose(f0);
  }

  return idx;
}

void cm_lookup( FLOAT3 v0, int cm_side, int *idx, float *val )
{
  FLOAT3 v = v0;

  int plane;
  float sx, sy;

  int major_plane =0;
  if( fabs(v.x) < fabs(v.y) )
    major_plane = 1;

  if( major_plane )
  {
    if( fabs(v.y) < fabs(v.z) )
      major_plane = 2;
  }else
  {
    if( fabs(v.x) < fabs(v.z) )
      major_plane = 2;
  }


  switch( major_plane )
  {
    case 0:
      if( v.x > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.x;
        sx = (1 + v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 0;
      }else
      {
        //{ -1,  1,  1 },  // v2
        v = v/(-v.x);
        sx = (1 - v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 1;
      }
      break;

    case 1:
      if( v.y > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.y;
        sx = (1 - v.x) / 2;
        sy = (1 + v.z) / 2;
        plane = 2;
      }else
      {
        //{  1, -1,  1 },  // v7
        v = v/(-v.y);
        sx = (1 - v.x) / 2;
        sy = (1 - v.z) / 2;
        plane = 3;
      }
      break;

    case 2:
      if( v.z > 0 )
      {
        //{  1,  1,  1 },  // v3
        v = v/v.z;
        sx = (1 - v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 4;
      }else
      {
        //{ -1,  1, -1 },  // v1
        v = v/(-v.z);
        sx = (1 + v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 5;
      }
      break;
  };


  sx = G_CLAMP( sx, 0, 1-FLT_EPSILON ) * cm_side + .5f;
  sy = G_CLAMP( sy, 0, 1-FLT_EPSILON ) * cm_side + .5f;

  int x0,y0, x1,y1;
    x0 = (int)G_CLAMP( floor(sx), 1, cm_side ) - 1;
    y0 = (int)G_CLAMP( floor(sy), 1, cm_side ) - 1;
    x1 = (int)G_CLAMP(  ceil(sx), 1, cm_side ) - 1;
    y1 = (int)G_CLAMP(  ceil(sy), 1, cm_side ) - 1;

  int idx0, idx1, idx2, idx3;
    idx0 = cm_side*cm_side*plane + y0*cm_side + x0;
    idx1 = cm_side*cm_side*plane + y0*cm_side + x1;
    idx2 = cm_side*cm_side*plane + y1*cm_side + x0;
    idx3 = cm_side*cm_side*plane + y1*cm_side + x1;

  float rx, ry;
    //sx = G_CLAMP( sx, 1, cm_side ) - 1;
    //sy = G_CLAMP( sy, 1, cm_side ) - 1;
    rx = sx - floorf(sx);
    ry = sy - floorf(sy);

  //FLOAT3 val;
  //  val = 
  //  ( (1-rx)*cm_data[idx0] + rx*cm_data[idx1] ) * (1-ry) +
  //  ( (1-rx)*cm_data[idx2] + rx*cm_data[idx3] ) * ry;
  //return val;

  idx[0] = idx0;
  idx[1] = idx1;
  idx[2] = idx2;
  idx[3] = idx3;
  val[0] = (1-rx)*(1-ry);
  val[1] = (  rx)*(1-ry);
  val[2] = (1-rx)*(  ry);
  val[3] = (  rx)*(  ry);
}

void em_lookup( FLOAT3 v0, int em_side, int *idx, float *val )
{
  FLOAT3 v = v0;

  int plane;
  float sx, sy;

  int major_plane =0;
  if( fabs(v.x) < fabs(v.y) )
    major_plane = 1;

  if( major_plane )
  {
    if( fabs(v.y) < fabs(v.z) )
      major_plane = 2;
  }else
  {
    if( fabs(v.x) < fabs(v.z) )
      major_plane = 2;
  }


  switch( major_plane )
  {
    case 0:
      if( v.x > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.x;
        sx = (1 + v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 0;
      }else
      {
        //{ -1,  1,  1 },  // v2
        v = v/(-v.x);
        sx = (1 - v.z) / 2;
        sy = (1 - v.y) / 2;
        plane = 1;
      }
      break;

    case 1:
      if( v.y > 0 )
      {
        //{  1,  1, -1 },  // v0
        v = v/v.y;
        sx = (1 - v.x) / 2;
        sy = (1 + v.z) / 2;
        plane = 2;
      }else
      {
        //{  1, -1,  1 },  // v7
        v = v/(-v.y);
        sx = (1 - v.x) / 2;
        sy = (1 - v.z) / 2;
        plane = 3;
      }
      break;

    case 2:
      if( v.z > 0 )
      {
        //{  1,  1,  1 },  // v3
        v = v/v.z;
        sx = (1 - v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 4;
      }else
      {
        //{ -1,  1, -1 },  // v1
        v = v/(-v.z);
        sx = (1 + v.x) / 2;
        sy = (1 - v.y) / 2;
        plane = 5;
      }
      break;
  };


  sx = G_CLAMP( sx, 0, 1 ) * (em_side-2) + .5f;
  sy = G_CLAMP( sy, 0, 1 ) * (em_side-2) + .5f;

  int x0,y0, x1,y1;
    x0 = (int)G_CLAMP( floor(sx), 0, em_side-1 );
    y0 = (int)G_CLAMP( floor(sy), 0, em_side-1 );
    x1 = (int)G_CLAMP(  ceil(sx), 0, em_side-1 );
    y1 = (int)G_CLAMP(  ceil(sy), 0, em_side-1 );

  int idx0, idx1, idx2, idx3;
    idx0 = em_side*em_side*plane + y0*em_side + x0;
    idx1 = em_side*em_side*plane + y0*em_side + x1;
    idx2 = em_side*em_side*plane + y1*em_side + x0;
    idx3 = em_side*em_side*plane + y1*em_side + x1;

  float rx, ry;
    rx = sx - floorf(sx);
    ry = sy - floorf(sy);

  idx[0] = idx0;
  idx[1] = idx1;
  idx[2] = idx2;
  idx[3] = idx3;
  val[0] = (1-rx)*(1-ry);
  val[1] = (  rx)*(1-ry);
  val[2] = (1-rx)*(  ry);
  val[3] = (  rx)*(  ry);
}

int sm_lookup( FLOAT3 v0, int smw, int smh )
{
  float theta, phi;
  Float32Angle( &theta, &phi, &v0 );
  if( phi<0 ) phi+= 2*G_PI;

  float sx, sy;
  sx = phi / ( 2 * G_PI );
  sy = theta / G_PI;

  sx = G_CLAMP( sx, 0, 1-FLT_EPSILON ) * smw  + .5f;
  sy = G_CLAMP( sy, 0, 1-FLT_EPSILON ) * smh + .5f;

  int x0,y0;
    x0 = (int)G_CLAMP( floor(sx), 1, smw  ) - 1;
    y0 = (int)G_CLAMP( floor(sy), 1, smh ) - 1;

  int idx0;
    idx0 = y0*smw + x0;

  return idx0;
}


void sm_lookup( FLOAT3 v0, int smw, int smh, int *idx, float *val )
{
  float theta, phi;
  Float32Angle( &theta, &phi, &v0 );
  if( phi<0 ) phi+= 2*G_PI;

  float sx, sy;
  sx = phi / ( 2 * G_PI );
  sy = theta / G_PI;

  sx = G_CLAMP( sx, 0, 1-FLT_EPSILON ) * smw  + .5f;
  sy = G_CLAMP( sy, 0, 1-FLT_EPSILON ) * smh + .5f;

  int x0,y0, x1,y1;
    x0 = (int)G_CLAMP( floor(sx), 1, smw  ) - 1;
    x1 = (int)G_CLAMP(  ceil(sx), 1, smw  ) - 1;
    y0 = (int)G_CLAMP( floor(sy), 1, smh ) - 1;
    y1 = (int)G_CLAMP(  ceil(sy), 1, smh ) - 1;

  int idx0, idx1, idx2, idx3;
    idx0 = y0*smw + x0;
    idx1 = y0*smw + x1;
    idx2 = y1*smw + x0;
    idx3 = y1*smw + x1;

  float rx, ry;
    rx = sx - floorf(sx);
    ry = sy - floorf(sy);

  idx[0] = idx0;
  idx[1] = idx1;
  idx[2] = idx2;
  idx[3] = idx3;
  val[0] = (1-rx)*(1-ry);
  val[1] = (  rx)*(1-ry);
  val[2] = (1-rx)*(  ry);
  val[3] = (  rx)*(  ry);
}

void get_hammersley( int n_point, GPfm &hsdir )
{
  int     k, kk;
  float  p, phi;

  float cost, sint;
  float cosp, sinp;

  hsdir.load( n_point, 1 );

  for( k=0; k<n_point; k++ )
  {
    cost = 0;
    for( kk = k, p = 0.5; kk; p *= 0.5, kk >>= 1 )
      if( kk & 1 )  // kk mod 2 == 1
        cost += p;

    // a slight shift and then map to [0, 2 pi)
    phi = 2*G_PI * ( (k+0.5f)/n_point );
    
    // map from [0,1] to [-1,1]
    cost = 2 * cost - 1;

    sint = sqrtf(  g_max( 1-cost*cost,0.f )  );
    sinp = sinf(phi);
    cosp = cosf(phi);
    
    hsdir.fm[k] = vnormalize( FLOAT3( sint*cosp, cost, -sint*sinp ) );
  }
}

void cm2dp( const GPfm &cmsrc, GPfm &dmsrc, int dpside )
{
  int cmside, dsside;
  float dc, ds;
  FLOAT3 v0, v1, v2;
  int i, j;
  float x, y;
  
  cmside = cmsrc.w;
  x = 2*float(cmside/2)/(cmside+1)-1;
  y = 2*float(cmside/2+1)/(cmside+1)-1;
  v0 = vnormalize( FLOAT3(x,x,1) );
  v1 = vnormalize( FLOAT3(y,x,1) );
  v2 = vnormalize( FLOAT3(x,y,1) );
  dc = area_spherical_triangle(v0,v1,v2);

  dsside = cmside;
  while(1)
  {
    x = 2*float(dsside/2)/(dsside+1)-1;
    y = 2*float(dsside/2+1)/(dsside+1)-1;
    v0 = vnormalize( FLOAT3(2*x,2*x,1-x*x-x*x) );
    v1 = vnormalize( FLOAT3(2*x,2*y,1-x*x-y*y) );
    v2 = vnormalize( FLOAT3(2*y,2*x,1-x*x-y*y) );
    ds = area_spherical_triangle(v0,v1,v2);
    if(10*ds<dc)
      break;
    dsside++;
  }

  GPfm emsrc, dpface;
  prepare_emsrc( cmsrc, emsrc );
  dmsrc.load( 2*dpside, dpside );

  dpface.load( dsside, dsside );
  for( j=0; j<dsside; j++ )
    for( i=0; i<dsside; i++ )
    {
      x = 1-2*(i+.5f)/dsside;
      y = 1-2*(j+.5f)/dsside;
      dpface.pm[j][i] =  em_lookup( vnormalize(FLOAT3(2*x,2*y,1-x*x-y*y)), emsrc.fm, emsrc.w );
    }
  dpface.resample( dpside, dpside );

  //for( j=0; j<dpface.h; j++ )
  //  for( i=0; i<dpface.w; i++ )
  //  {
  //    x = 1-2*(i+.5f)/dpface.w;
  //    y = 1-2*(j+.5f)/dpface.h;
  //    if( x*x+y*y>1 )
  //    dpface.pm[j][i] = 0;
  //  }

  dmsrc.draw( dpface, 0,0 );

  dpface.load( dsside, dsside );
  for( j=0; j<dsside; j++ )
    for( i=0; i<dsside; i++ )
    {
      x = 1-2*(i+.5f)/dsside;
      y = 1-2*(j+.5f)/dsside;
      dpface.pm[j][i] =  em_lookup( vnormalize(FLOAT3(-2*x,2*y,-(1-x*x-y*y))), emsrc.fm, emsrc.w );
    }

  dpface.resample( dpside, dpside );

  //for( j=0; j<dpface.h; j++ )
  //  for( i=0; i<dpface.w; i++ )
  //  {
  //    x = 1-2*(i+.5f)/dpface.w;
  //    y = 1-2*(j+.5f)/dpface.h;
  //    if( x*x+y*y>1 )
  //    dpface.pm[j][i] = 0;
  //  }
  dmsrc.draw( dpface, dpside,0 );
}


