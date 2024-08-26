#ifndef G_OBJ_H
#define G_OBJ_H

//#include <GL/glew.h>
//#include <GL/gl.h>

#include "g_common.h"
#include "g_pfm.h"

#define GOBJ_VERTEX 1
#define GOBJ_TEXCOR 1<<1
#define GOBJ_NORMAL 1<<2
#define GOBJ_TANGENT 1<<3
#define GOBJ_NO_MATERIAL 1<<4

#define GOBJ_V__ (GOBJ_VERTEX)
#define GOBJ_VT_ (GOBJ_VERTEX | GOBJ_TEXCOR)
#define GOBJ_V_N (GOBJ_VERTEX | GOBJ_NORMAL)
#define GOBJ_VTN (GOBJ_VERTEX | GOBJ_TEXCOR | GOBJ_NORMAL)
#define GOBJ_VTNT (GOBJ_VERTEX | GOBJ_TEXCOR | GOBJ_NORMAL | GOBJ_TANGENT)

typedef struct _GEdge
{
  int n0, v0, t0;
  int n1, v1, t1;
}GEdge;

class Group
{
  public:
    char name[256];
    int n_face;
    int type;

    int n_buf;
    int *tri;
    void add( int face_idx );
    void remove( int face_idx );
    void clear();

    Group();
    ~Group();
};

class GObj
{
  public:
    GObj();
    ~GObj();
    void obj_merge( const GObj &a, const GObj &b );
    void obj_merge( const GObj &b );  // depreciated
    void merge( const GObj &b );
    FLOAT4 unitize();
    void load_vertex( int number_of_vertex );
    void load_normal( int number_of_normal );
    void load_texc( int number_of_texc );
    void load_face( int number_of_face );
    void load_mgroup( int number_of_group );
    void load_sgroup( int number_of_group );
    void load_mtl( const char *spath=0, int tex_size=512 );
    void load( const char *spath=0, int tex_size=512 );
    void save( const char *spath ) const;
    void at_least_one_mgroup();

    void texfilter( unsigned int filter_type );
    void zero_center();

    void clear();

    void draw( int draw_option = GOBJ_VTN ) const;
    unsigned int list( int draw_option = GOBJ_VTN );
    unsigned int list_by_group( int draw_option = GOBJ_VTN );
    void calculate_vertex_normal( float angle );
    void calculate_face_normal();
    double calculate_face_area( int ti );
    void calculate_tangent_space();
    void translate( float x, float y, float z );
    void rotate_euler( float rx, float ry, float rz );
    void rotate( GQuat q );
    void rotate( float angle, float x, float y, float z );
    void scale( float scale );
    void reversewinding();
    FLOAT3 get_face_normal( int face_index ) const;
    FLOAT3 get_face_normal( FLOAT3 v0, FLOAT3 v1, FLOAT3 v2 ) const;

    void unitize_normal();
    void clean_redundant_faces();

    void trim();

    void weld();
    void weld_normal();

    typedef struct _Material
    {
      char name[256];           // name of material
      FLOAT4 ambient;           // ambient component
      FLOAT4 diffuse;           // diffuse component
      FLOAT4 specular;          // specular component
      float shininess;          // specular exponent
      FLOAT4 ward;
      float bump;
      int illum;
      unsigned int texid;
      char texname[256];
    }Material;

    int n_material;
    Material *material;

    Group default_group;

    int n_mgroup;
    Group *mgroup;

    int n_sgroup;
    Group *sgroup;

    int n_vertex;
    int n_normal;
    int n_texc;

    GPfm vertex;
    GPfm normal;
    GPfm texc;
    GPfm tangent;
    GPfm binormal;

    int n_face;
    int *face_nidx;
    int *face_vidx;
    int *face_tidx;
    int *face_eidx;

    unsigned int tex0;

    bool favour_sgroup;
    bool modified;


    char obj_path[256];
    char mtl_path[256];
    char tex_path[256];
    void printf_info();



    GEdge *s_edge;
    int n_edge;
    bool closed;
    void calculate_edge_info();
    void save_edge_info( const char *spath );
    void load_edge_info( const char *spath );

    float calculate_smoothness();
    void force_single_sided();
    void calculate_access( int *&a0, int **&b0, const int *access, int n_access );
    Material getmtl( const char *id ) const;
    void setmtl( const char *id, const Material &mtl, int tex_size=512 );
    void usemtl( const Material &mtl ) const;
    void save_mtl( const char *spath ) const;

    void unfold_in_normal_access();
    void draw_point( int vertex_index );
    void draw_point_set( const GStack<int> &point_set );
    void draw_edge( int eidx );
    void draw_face( int tidx );
    void draw_group( const Group &grp, int draw_option=GOBJ_VTN ) const;
    void gen_vertex_array( const Group &grp, int &ndat, float* &dat ) const;

  private:
    void loadx( const char *spath );
    void calculate_face_tangent( GPfm &face_tangent );
    void openwrl( const char *spath, int tex_size=512 );
};

#endif