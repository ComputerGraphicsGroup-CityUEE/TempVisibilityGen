#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <float.h>
#include <math.h>
#include "g_vector.h"


DOUBLE3 DOUBLE3::operator+(const DOUBLE3 &a) const { return DOUBLE3(x + a.x, y + a.y, z + a.z); }
DOUBLE3 DOUBLE3::operator-(const DOUBLE3 &a) const { return DOUBLE3(x - a.x, y - a.y, z - a.z); }
DOUBLE3 DOUBLE3::operator*(const DOUBLE3 &a) const { return DOUBLE3(x * a.x, y * a.y, z * a.z); }
DOUBLE3 DOUBLE3::operator/(const DOUBLE3 &a) const { return DOUBLE3(x / a.x, y / a.y, z / a.z); }

DOUBLE3 DOUBLE3::operator+(const double &a) const { return DOUBLE3(x + a, y + a, z + a); }
DOUBLE3 DOUBLE3::operator-(const double &a) const { return DOUBLE3(x - a, y - a, z - a); }
DOUBLE3 DOUBLE3::operator*(const double &a) const { return DOUBLE3(x * a, y * a, z * a); }
DOUBLE3 DOUBLE3::operator/(const double &a) const { return DOUBLE3(x / a, y / a, z / a); }

DOUBLE3& DOUBLE3::operator+=(const DOUBLE3 &a)
{
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}
DOUBLE3& DOUBLE3::operator-=(const DOUBLE3 &a)
{
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}
DOUBLE3& DOUBLE3::operator*=(const DOUBLE3 &a)
{
  x *= a.x;
  y *= a.y;
  z *= a.z;
  return *this;
}
DOUBLE3& DOUBLE3::operator/=(const DOUBLE3 &a)
{
  x /= a.x;
  y /= a.y;
  z /= a.z;
  return *this;
}

DOUBLE3 operator+(double a, const DOUBLE3 &b) { return DOUBLE3(a + b.x, a + b.y, a + b.z); }
DOUBLE3 operator-(double a, const DOUBLE3 &b) { return DOUBLE3(a - b.x, a - b.y, a - b.z); }
DOUBLE3 operator*(double a, const DOUBLE3 &b) { return DOUBLE3(a * b.x, a * b.y, a * b.z); }
DOUBLE3 operator/(double a, const DOUBLE3 &b) { return DOUBLE3(a / b.x, a / b.y, a / b.z); }

DOUBLE3 DOUBLE3::operator-() const { return DOUBLE3(-x, -y, -z); }
DOUBLE3 DOUBLE3::operator+() const { return DOUBLE3(x, y, z); }

DOUBLE3& DOUBLE3::operator=(double a) { x = a; y = a; z = a; return *this; }
DOUBLE3& DOUBLE3::operator=(const FLOAT3 &a) { x = a.x; y = a.y; z = a.z; return *this; }

// row vector affine transform, 
DOUBLE3 DOUBLE3::rmul(const double *m) const
{
  double tx, ty, tz, tw;
  tx = m[0] * x + m[4] * y + m[8] * z + m[12];
  ty = m[1] * x + m[5] * y + m[9] * z + m[13];
  tz = m[2] * x + m[6] * y + m[10] * z + m[14];
  tw = m[3] * x + m[7] * y + m[11] * z + m[15];
  return DOUBLE3(tx / tw, ty / tw, tz / tw);
}
double DOUBLE3::norm() const
{
  return sqrt(x*x + y*y + z*z);
}

bool DOUBLE3::operator==(const DOUBLE3 &b) const
{
  return x == b.x && y == b.y && z == b.z;
}
bool DOUBLE3::operator!=(const DOUBLE3 &b) const
{
  return x != b.x || y != b.y || z != b.z;
}

double vdot(const DOUBLE3 &a, const DOUBLE3 &b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

DOUBLE3 vnormalize(const DOUBLE3 &a)
{
  double l = sqrt(vdot(a, a));
  if (l == 0)
    return 0;
  else
    return a / l;
}

DOUBLE3 vcross(const DOUBLE3 &a, const DOUBLE3 &b)
{
  DOUBLE3 c;
  c.x = a.y*b.z - b.y*a.z;
  c.y = -a.x*b.z + b.x*a.z;
  c.z = a.x*b.y - b.x*a.y;
  return c;
}

void Angle2Double3(const double theta, const double phi, DOUBLE3 *v)	//convert a solid angle to a vector
{
  double cost, sint, cosp, sinp;
  cost = cos(theta);
  sint = sin(theta);
  cosp = cos(phi);
  sinp = sin(phi);

  v->x = sint * cosp;
  v->y = cost;
  v->z = -sint * sinp;
}

void Double32Angle(double *theta, double *phi, const DOUBLE3 *v)
{
  double r = sqrt(vdot(*v, *v)) + DBL_EPSILON;
  *theta = acos(v->y / r);
  *phi = atan2(-v->z / r, v->x / r + DBL_EPSILON);
}

double fmax( const DOUBLE3 &a )
{
  double res = a.x;
  res = res > a.y ? res : a.y; 
  res = res > a.z ? res : a.z; 
  return res;
}


FLOAT3 FLOAT3::operator+( const FLOAT3 &a ) const { return FLOAT3( x+a.x, y+a.y, z+a.z ); }
FLOAT3 FLOAT3::operator-( const FLOAT3 &a ) const { return FLOAT3( x-a.x, y-a.y, z-a.z ); }
FLOAT3 FLOAT3::operator*( const FLOAT3 &a ) const { return FLOAT3( x*a.x, y*a.y, z*a.z ); }
FLOAT3 FLOAT3::operator/( const FLOAT3 &a ) const { return FLOAT3( x/a.x, y/a.y, z/a.z ); }

FLOAT3 FLOAT3::operator+( const float &f ) const { return FLOAT3( x+f, y+f, z+f ); }
FLOAT3 FLOAT3::operator-( const float &f ) const { return FLOAT3( x-f, y-f, z-f ); }
FLOAT3 FLOAT3::operator*( const float &f ) const { return FLOAT3( x*f, y*f, z*f ); }
FLOAT3 FLOAT3::operator/( const float &f ) const { return FLOAT3( x/f, y/f, z/f ); }

FLOAT3& FLOAT3::operator+=( const FLOAT3 &a ){ this->x+=a.x; this->y+=a.y; this->z+=a.z; return *this; }
FLOAT3& FLOAT3::operator-=( const FLOAT3 &a ){ this->x-=a.x; this->y-=a.y; this->z-=a.z; return *this; }
FLOAT3& FLOAT3::operator*=( const FLOAT3 &a ){ this->x*=a.x; this->y*=a.y; this->z*=a.z; return *this; }
FLOAT3& FLOAT3::operator/=( const FLOAT3 &a ){ this->x/=a.x; this->y/=a.y; this->z/=a.z; return *this; }

FLOAT3 operator*( float a, const FLOAT3 &b ){ return FLOAT3( a*b.x, a*b.y, a*b.z ); }
FLOAT3 operator+( float a, const FLOAT3 &b ){ return FLOAT3( a+b.x, a+b.y, a+b.z ); }
FLOAT3 operator-( float a, const FLOAT3 &b ){ return FLOAT3( a-b.x, a-b.y, a-b.z ); }
FLOAT3 operator/( float a, const FLOAT3 &b ){ return FLOAT3( a/b.x, a/b.y, a/b.z ); }


FLOAT3 FLOAT3::operator-() const { return FLOAT3( -x,-y,-z ); }


FLOAT3::FLOAT3( const BYTE3 &a )
{
  this->x = (float)a.x;
  this->y = (float)a.y;
  this->z = (float)a.z;
}

FLOAT3& FLOAT3::operator=( const BYTE3 &a )
{
  this->x = (float)a.x;
  this->y = (float)a.y;
  this->z = (float)a.z;

  return *this;
}

FLOAT3& FLOAT3::operator=( float a )
{
  this->x = a;
  this->y = a;
  this->z = a;

  return *this;
}

FLOAT3 FLOAT3::rmul( const float *m ) const		//vector right multiplication, m is a 4x4 matrix, column major  
{
  FLOAT3 v;
    v = FLOAT3(
      x * m[0] + y * m[4] + z * m[8], 
      x * m[1] + y * m[5] + z * m[9], 
      x * m[2] + y * m[6] + z * m[10]
      );
  return v;
}

float FLOAT3::norm() const
{
  return sqrtf( vdot( this,this ) );
}

bool FLOAT3::operator==( const FLOAT3 &b ) const
{
  if( x==b.x && y==b.y && z==b.z )
    return true;
  return false;
}

bool FLOAT3::operator!=( const FLOAT3& b ) const
{
  if( x!=b.x || y!=b.y || z!=b.z )
    return true;
  return false;
}


FLOAT3 f3abs( const FLOAT3 &a )
{
  FLOAT3 c; 
    c.x = fabsf(a.x);
    c.y = fabsf(a.y);
    c.z = fabsf(a.z);
  return c;
}

float vangle( const FLOAT3 &a, const FLOAT3 &b )
{
  float d;
  d = vdot(a,b);
  d = d<-1?-1:(d<1?d:1);
  return acosf(d);
}

FLOAT4 FLOAT4::operator+( const FLOAT4 &a ) const { return FLOAT4( x+a.x, y+a.y, z+a.z, w+a.w ); }
FLOAT4 FLOAT4::operator-( const FLOAT4 &a ) const { return FLOAT4( x-a.x, y-a.y, z-a.z, w-a.w ); }
FLOAT4 FLOAT4::operator*( const FLOAT4 &a ) const { return FLOAT4( x*a.x, y*a.y, z*a.z, w*a.w ); }
FLOAT4 FLOAT4::operator/( const FLOAT4 &a ) const { return FLOAT4( x/a.x, y/a.y, z/a.z, w/a.w ); }

FLOAT4 FLOAT4::operator+( const float &f ) const { return FLOAT4( x+f, y+f, z+f, w+f ); }
FLOAT4 FLOAT4::operator-( const float &f ) const { return FLOAT4( x-f, y-f, z-f, w-f ); }
FLOAT4 FLOAT4::operator*( const float &f ) const { return FLOAT4( x*f, y*f, z*f, w*f ); }
FLOAT4 FLOAT4::operator/( const float &f ) const { return FLOAT4( x/f, y/f, z/f, w/f ); }

FLOAT4& FLOAT4::operator+=( const FLOAT4 &a ){ this->x+=a.x; this->y+=a.y; this->z+=a.z; this->w+=a.w; return *this; }
FLOAT4& FLOAT4::operator-=( const FLOAT4 &a ){ this->x-=a.x; this->y-=a.y; this->z-=a.z; this->w-=a.w; return *this; }
FLOAT4& FLOAT4::operator*=( const FLOAT4 &a ){ this->x*=a.x; this->y*=a.y; this->z*=a.z; this->w*=a.w; return *this; }
FLOAT4& FLOAT4::operator/=( const FLOAT4 &a ){ this->x/=a.x; this->y/=a.y; this->z/=a.z; this->w/=a.w; return *this; }

FLOAT4 operator+( float a, const FLOAT4 &b ){ return FLOAT4( a+b.x, a+b.y, a+b.z, a+b.w ); }
FLOAT4 operator-( float a, const FLOAT4 &b ){ return FLOAT4( a-b.x, a-b.y, a-b.z, a-b.w ); }
FLOAT4 operator*( float a, const FLOAT4 &b ){ return FLOAT4( a*b.x, a*b.y, a*b.z, a*b.w ); }
FLOAT4 operator/( float a, const FLOAT4 &b ){ return FLOAT4( a/b.x, a/b.y, a/b.z, a/b.w ); }

FLOAT4 FLOAT4::operator-() const { return FLOAT4( -x,-y,-z,-w ); }

bool FLOAT4::operator==( const FLOAT4 &b ) const
{
  if( x==b.x && y==b.y && z==b.z && w==b.w )
    return true;
  return false;
}
bool FLOAT4::operator!=( const FLOAT4 &b ) const
{
  if( x!=b.x || y!=b.y || z!=b.z || w!=b.w )
    return true;
  return false;
}

FLOAT4 FLOAT4::operator*( const float *mx ) const
{
  FLOAT4 hv;
  hv.x = vdot( (float*)this, &mx[0], 4 );
  hv.y = vdot( (float*)this, &mx[4], 4 );
  hv.z = vdot( (float*)this, &mx[8], 4 );
  hv.w = vdot( (float*)this, &mx[12], 4 );
  return hv;
}
FLOAT4 operator*( const float *mx, const FLOAT4 &b )
{
  FLOAT4 hv;
  hv.x = mx[0]*b.x + mx[4]*b.y + mx[ 8]*b.z + mx[12]*b.w;
  hv.y = mx[1]*b.x + mx[5]*b.y + mx[ 9]*b.z + mx[13]*b.w;
  hv.z = mx[2]*b.x + mx[6]*b.y + mx[10]*b.z + mx[14]*b.w;
  hv.w = mx[3]*b.x + mx[7]*b.y + mx[11]*b.z + mx[15]*b.w;
  return hv;
}

FLOAT4 FLOAT4::rmul( const float *m ) const
{
  FLOAT4 v;
    v = FLOAT4(
      x*m[0] + y*m[4] + z*m[ 8] + w*m[12], 
      x*m[1] + y*m[5] + z*m[ 9] + w*m[13], 
      x*m[2] + y*m[6] + z*m[10] + w*m[14],
      x*m[3] + y*m[7] + z*m[11] + w*m[15]
      );
  return v;
}

FLOAT4::FLOAT4( const BYTE4 &a )
{
  this->x = (float)a.x;
  this->y = (float)a.y;
  this->z = (float)a.z;
  this->w = (float)a.w;
}

FLOAT4& FLOAT4::operator=( const BYTE4 &a )
{
  this->x = (float)a.x;
  this->y = (float)a.y;
  this->z = (float)a.z;
  this->w = (float)a.w;

  return *this;
}






BYTE3& BYTE3::operator=( const GBYTE &a )
{
  this->x = a;
  this->y = a;
  this->z = a;
  return *this;
}

#define G_CLAMP(x,a,b)( ((x)<(a))? (a) : (((x)<(b))?(x):(b))  )
BYTE3& BYTE3::operator=( const FLOAT3 &a )
{
  this->x = (GBYTE)G_CLAMP( a.x, 0, 255 );
  this->y = (GBYTE)G_CLAMP( a.y, 0, 255 );
  this->z = (GBYTE)G_CLAMP( a.z, 0, 255 );
  return *this;
}

bool BYTE3::operator==( const BYTE3 &b ) const
{
  if( x==b.x && y==b.y && z==b.z )
    return true;
  return false;
}


float vdot(const FLOAT3 &a, const FLOAT3 &b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

float vdot( const float *a, const float *b, int n )
{
  int i;
  float c;
  for( i=0, c=0; i<n; i++, a++, b++ )
    c += *a * *b;
  return c;
}


float vdot(const FLOAT3 *a, const FLOAT3 *b)
{
  return a->x*b->x + a->y*b->y + a->z*b->z;
}

FLOAT3 vmin(const FLOAT3 &a, const FLOAT3 &b)
{
  FLOAT3 c;;
    c.x = (a.x < b.x) ? a.x : b.x;
    c.y = (a.y < b.y) ? a.y : b.y;
    c.z = (a.z < b.z) ? a.z : b.z;
  return c;
}

FLOAT4 vmin(const FLOAT4 &a, const FLOAT4 &b)
{
  FLOAT4 c;;
    c.x = (a.x < b.x) ? a.x : b.x;
    c.y = (a.y < b.y) ? a.y : b.y;
    c.z = (a.z < b.z) ? a.z : b.z;
    c.w = (a.w < b.w) ? a.w : b.w;
  return c;
}

FLOAT3 vmax(const FLOAT3 &a, const FLOAT3 &b)
{
  FLOAT3 c;;
    c.x = (a.x > b.x) ? a.x : b.x;
    c.y = (a.y > b.y) ? a.y : b.y;
    c.z = (a.z > b.z) ? a.z : b.z;
  return c;
}

FLOAT4 vmax(const FLOAT4 &a, const FLOAT4 &b)
{
  FLOAT4 c;;
    c.x = (a.x > b.x) ? a.x : b.x;
    c.y = (a.y > b.y) ? a.y : b.y;
    c.z = (a.z > b.z) ? a.z : b.z;
    c.w = (a.w > b.w) ? a.w : b.w;
  return c;
}

FLOAT3 vcross( const FLOAT3 &a, const FLOAT3 &b )
{
  FLOAT3 c;
    c.x =  a.y*b.z - b.y*a.z;
    c.y = -a.x*b.z + b.x*a.z;
    c.z =  a.x*b.y - b.x*a.y;

  return c;
}

// Halve arc between unit vectors v0 and v1.
FLOAT3 vbisect( const FLOAT3 &a, const FLOAT3 &b )
{
  FLOAT3 v = a+b;
  float Nv = vdot(&v,&v);

  if (Nv < 1.0e-5) 
    v = FLOAT3(0, 0, 1);
  else
    v = v/sqrtf(Nv);

  return v;
}


void vnormalize( FLOAT3 *a )
{
  float l = sqrtf( vdot(a,a) );

  if( l==0 )
    *a = 0;
  else
    *a = *a/l;
}

FLOAT3 vnormalize( const FLOAT3 &a )
{
  float l = sqrtf( vdot(a,a) );
  if( l==0 )
    return 0;
  else
    return a/l;
}

float vnormalize( float *a, int n )
{
  int j;
  float d;
  d = sqrtf(vdot(a,a,n));
  for( j=0; j<n; j++ )
    a[j] = a[j]/d;
  return d;
}


FLOAT3 intersect( const FLOAT3 &u0, const FLOAT3 &u1, const FLOAT3 &v0, const FLOAT3 &v1 )
{
  FLOAT3 du, dv;
    du = u1 - u0;
    dv = v1 - v0;
    return vcross( v0-u0, dv ).norm() / vcross( du,dv ).norm() * du + u0;
}

FLOAT3 proj_onto_plane( const FLOAT3 &v, const FLOAT3 &n, const FLOAT3 &v0 )
{
  return v - vdot( v - v0, n ) * n;
}

FLOAT3 sqrtf3( const FLOAT3 &a )
{
  return FLOAT3( sqrtf(a.x), sqrtf(a.y), sqrtf(a.z) );
}

void Angle2Float3( const float theta, const float phi, FLOAT3 *v)	//convert a solid angle to a vector
{
  float cost, sint, cosp, sinp;
    cost  = cosf(theta);
    sint  = sinf(theta);
    cosp  = cosf(phi);
    sinp  = sinf(phi);
  
  v->x = sint * cosp;
  v->y = cost;
  v->z = -sint * sinp;
}

void Float32Angle( float *theta, float *phi, const FLOAT3 *v)
{
  float r = sqrtf( vdot(v,v) )+FLT_EPSILON;
  *theta  = acosf( v->y/r );
  *phi    = atan2f( -v->z/r, v->x/r+FLT_EPSILON );
}

void Mul_M33xC3( const float *m, const FLOAT3 *tX, FLOAT3 *tY )
{
  tY->x = vdot( (FLOAT3*)&m[0], tX );
  tY->y = vdot( (FLOAT3*)&m[3], tX );
  tY->z = vdot( (FLOAT3*)&m[6], tX );
}

void Mul_R3xM33( const float *m, const FLOAT3 *tX, FLOAT3 *tY )
{
  tY->x = tX->x * m[0] + tX->y * m[3] + tX->z * m[6] ;
  tY->y = tX->x * m[1] + tX->y * m[4] + tX->z * m[7] ;
  tY->z = tX->x * m[2] + tX->y * m[5] + tX->z * m[8] ;
}

void vsub( const float *a, const float *b, float *c, int vsize )
{
  int i;
  for( i=0; i<vsize; i++ )
    c[i] = a[i] - b[i];
}

float vnorm( const float *a, int vsize )
{
  float sum;
  int i;

  sum = 0;
  for( i=0; i<vsize; i++ )
    sum += a[i]*a[i];

  return sqrtf(sum);
}

void vperturb( const float *a, float *b, int vsize, float perturb_factor )
{
  int i;
  float al = vnorm(a,vsize);

  float *c = (float*) malloc( vsize*sizeof(float) );

  for( i=0; i<vsize; i++ )
    c[i] = 2*float( rand() )/ RAND_MAX - 1;

  float vecl;
    vecl = vnorm( c, vsize );

  for( i=0; i<vsize; i++ )
    b[i] = a[i] + c[i] / vecl * al * perturb_factor;

  free(c);
}





GQuat GQuat::operator~() const
{
  GQuat a;
    a.x = -x;
    a.y = -y;
    a.z = -z;
    a.w = w;
  return a;
}


GQuat::GQuat( const FLOAT3 &axis, const float &angle )
{
  x = sinf(angle/2) * axis.x;
  y = sinf(angle/2) * axis.y;
  z = sinf(angle/2) * axis.z;
  w = cosf(angle/2);
}

GQuat::GQuat( const FLOAT3 &v0, const FLOAT3 &v1 )
{
//  x = v0.y*v1.z - v0.z*v1.y;
//  y = v0.z*v1.x - v0.x*v1.z;
//  z = v0.x*v1.y - v0.y*v1.x;
//  w = v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
  FLOAT3 axis; float area;
  float angle; float len;
    axis = vcross( v0, v1 );
    area = axis.norm();
    len  = vdot( v0, v1 );
    len  = G_CLAMP( len, -1, 1 );
    angle = acosf(len);

  if( area==0 && len<0 )
  {
    if( fabsf(v0.x)<fabsf(v0.y) && fabsf(v0.x)<fabsf(v0.z) )
      axis = vcross( v0, FLOAT3(1,0,0) );
    else if( fabsf(v0.y)<fabsf(v0.z) )
      axis = vcross( v0, FLOAT3(0,1,0) );
    else
      axis = vcross( v0, FLOAT3(0,0,1) );
    angle = acosf(-1);
  }

  *this = GQuat( vnormalize(axis), angle );
}

GQuat GQuat::operator*( const GQuat& b ) const		//Hamilton product
{
  GQuat c;

    c.w = w*b.w - x*b.x - y*b.y - z*b.z;
    c.x = w*b.x + x*b.w + y*b.z - z*b.y;
    c.y = w*b.y - x*b.z + y*b.w + z*b.x;
    c.z = w*b.z + x*b.y - y*b.x + z*b.w;

  return c;
}

GQuat GQuat::operator+( const GQuat& b ) const
{
  GQuat c;

    c.x = x + b.x;
    c.y = y + b.y;
    c.z = z + b.z;
    c.w = w + b.w;

  return c;
}

GQuat GQuat::operator*( const float& b ) const
{
  GQuat c;
    c.x = x*b;
    c.y = y*b;
    c.z = z*b;
    c.w = w*b;

  return c;
}

GQuat GQuat::operator/( const float& b ) const
{
  GQuat c;
    c.x = x/b;
    c.y = y/b;
    c.z = z/b;
    c.w = w/b;

  return c;
}

GQuat GQuat::operator/( const GQuat& b ) const
{
  float n2 = b.x*b.x + b.y*b.y + b.z*b.z + b.w*b.w;
  return *this  *  ~b/n2;
}

float GQuat::norm() const
{
  return sqrtf( x*x + y*y + z*z + w*w );
}

GQuat GQuat::normalize() const
{
  float d;
  GQuat c;
    d = this->norm();
    if(d>0)
      c = *this / d;
  return c;
}


void GQuat::matrix( float *m ) const
{

  GQuat n = this->normalize();
  float xx, yy, zz;
  float xy, yz, zw;
  float xz, yw;
  float xw;

    xx = n.x*n.x;  yy = n.y*n.y;  zz = n.z*n.z;
    xy = n.x*n.y;  yz = n.y*n.z;  zw = n.z*n.w;
    xz = n.x*n.z;  yw = n.y*n.w;
    xw = n.x*n.w;
        
    memset( m, 0, 16*sizeof(float) );

    m[0]  = 1 - 2 * ( yy + zz );
    m[1]  =     2 * ( xy - zw );
    m[2]  =     2 * ( xz + yw );
  
    m[4]  =     2 * ( xy + zw );
    m[5]  = 1 - 2 * ( xx + zz );
    m[6]  =     2 * ( yz - xw );
  
    m[8]  =     2 * ( xz - yw );
    m[9]  =     2 * ( yz + xw );
    m[10] = 1 - 2 * ( xx + yy );

    m[15] = 1;
}

void GQuat::axis_angle( FLOAT3 &axis, float &angle ) const
{
  GQuat qn = this->normalize();

  float cos_t, sin_t;
    cos_t = qn.w;
    sin_t = sqrtf( 1 - cos_t * cos_t );
    angle = 2 * acosf( cos_t );

    if ( fabs( sin_t ) < 0.0005 )
    {
      axis.x = qn.x;
      axis.y = qn.y;
      axis.z = qn.z;
    }else
    {
      axis.x = qn.x/sin_t;
      axis.y = qn.y/sin_t;
      axis.z = qn.z/sin_t;
    }
}

bool GQuat::operator==( const GQuat& b ) const
{
  if( x == b.x &&  y == b.y &&  z == b.z &&  w == b.w )
    return true;
  return false;
}


///////////////////////////////////////////////////////////////////////////////////
//
// Quaternion Interpolation from
// http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/
//
GQuat GQuat::slerp( const GQuat &a, const GQuat &b, float t )
{
  GQuat qa, qb, a_b;
  float sin_t, cos_t, theta;
  float sa, sb;

    qa = a.normalize();
    qb = b.normalize();

    cos_t = qa.w*qb.w + qa.x*qb.x + qa.y*qb.y + qa.z*qb.z;
    cos_t = cos_t< 1 ? cos_t :  1 ;
    cos_t = cos_t>-1 ? cos_t : -1 ;
    sin_t = sqrtf(1-cos_t*cos_t);

    cos_t = G_CLAMP( cos_t, -1,1 );

	//if theta = 0, a_b can be a, b, or any linearly combined vector between a and b
	//if theta = 180, a_b can be any vector which is normal to a or b
	//these two condition can be joint to one condition: a and b are in the same line.
    if ( fabsf( sin_t ) < 0.0005 )
    {
      sa = 1-t;
      sb = t;
    }else
    {
      theta = acosf( cos_t );
      sa = sinf( theta * (1-t) )/sin_t;
      sb = sinf( theta * t )/sin_t;
    }

    a_b = qa * sa + qb * sb;

    return a_b;
}



FLOAT3 MouseOnSphere( int x0, int y0, int width, int height )
{
  float x, y, z, w;

    x = 2.f*(x0+.5f)/width-1;
    y = 2.f*((height-y0-1)+.5f)/height-1;
    w = x*x + y*y;

  if( w>1.0 ) 
  {
    float scale = 1.0f/sqrtf(w);
    x *= scale; 
    y *= scale;
    z = 0.0f;
  }else
  {
    z = sqrtf(1 - w);
  }

  return FLOAT3(x,y,z);
}


DOUBLE4 DOUBLE4::operator+(const DOUBLE4 &a) const { return DOUBLE4(x + a.x, y + a.y, z + a.z, w + a.w); }
DOUBLE4 DOUBLE4::operator-(const DOUBLE4 &a) const { return DOUBLE4(x - a.x, y - a.y, z - a.z, w - a.w); }
DOUBLE4 DOUBLE4::operator*(const DOUBLE4 &a) const { return DOUBLE4(x * a.x, y * a.y, z * a.z, w * a.w); }
DOUBLE4 DOUBLE4::operator/(const DOUBLE4 &a) const { return DOUBLE4(x / a.x, y / a.y, z / a.z, w / a.w); }

DOUBLE4 DOUBLE4::operator+(const double &a) const { return DOUBLE4(x + a, y + a, z + a, w + a); }
DOUBLE4 DOUBLE4::operator-(const double &a) const { return DOUBLE4(x - a, y - a, z - a, w - a); }
DOUBLE4 DOUBLE4::operator*(const double &a) const { return DOUBLE4(x * a, y * a, z * a, w * a); }
DOUBLE4 DOUBLE4::operator/(const double &a) const { return DOUBLE4(x / a, y / a, z / a, w / a); }

DOUBLE4& DOUBLE4::operator+=(const DOUBLE4 &a)
{
  x += a.x;
  y += a.y;
  z += a.z;
  w += a.w;
  return *this;
}
DOUBLE4& DOUBLE4::operator-=(const DOUBLE4 &a)
{
  x -= a.x;
  y -= a.y;
  z -= a.z;
  w -= a.w;
  return *this;
}
DOUBLE4& DOUBLE4::operator*=(const DOUBLE4 &a)
{
  x *= a.x;
  y *= a.y;
  z *= a.z;
  w *= a.w;
  return *this;
}
DOUBLE4& DOUBLE4::operator/=(const DOUBLE4 &a)
{
  x /= a.x;
  y /= a.y;
  z /= a.z;
  w /= a.w;
  return *this;
}

DOUBLE4 operator+(double a, const DOUBLE4 &b) { return DOUBLE4(a + b.x, a + b.y, a + b.z, a + b.w); }
DOUBLE4 operator-(double a, const DOUBLE4 &b) { return DOUBLE4(a - b.x, a - b.y, a - b.z, a - b.w); }
DOUBLE4 operator*(double a, const DOUBLE4 &b) { return DOUBLE4(a * b.x, a * b.y, a * b.z, a * b.w); }
DOUBLE4 operator/(double a, const DOUBLE4 &b) { return DOUBLE4(a / b.x, a / b.y, a / b.z, a / b.w); }

DOUBLE4 operator*(double* m, const DOUBLE4 &b)
{
  return DOUBLE4(
    m[0]*b.x + m[4]*b.y + m[ 8]*b.z + m[12]*b.w,
    m[1]*b.x + m[5]*b.y + m[ 9]*b.z + m[13]*b.w,
    m[2]*b.x + m[6]*b.y + m[10]*b.z + m[14]*b.w,
    m[3]*b.x + m[7]*b.y + m[11]*b.z + m[15]*b.w );
}

DOUBLE4 DOUBLE4::operator-() const { return DOUBLE4(-x, -y, -z, -w); }
DOUBLE4 DOUBLE4::operator+() const { return *this; }

DOUBLE4& DOUBLE4::operator=(double a) { x = a; y = a; z = a; w = a; return *this; }

double DOUBLE4::norm() const
{
  return sqrt(x*x + y*y + z*z + w*w);
}

bool DOUBLE4::operator==(const DOUBLE4 &b) const
{
  return x == b.x && y == b.y && z == b.z && w == b.w;
}
bool DOUBLE4::operator!=(const DOUBLE4 &b) const
{
  return x != b.x || y != b.y || z != b.z || w != b.w;
}

DOUBLE4 log( const DOUBLE4 &a ){ return DOUBLE4( log(a.x), log(a.y), log(a.z), log(a.w) ); }
DOUBLE4 tanh( const DOUBLE4 &a ){ return DOUBLE4( tanh(a.x), tanh(a.y), tanh(a.z), tanh(a.w) ); }
DOUBLE4 fabs( const DOUBLE4 &a ){ return DOUBLE4( fabs(a.x), fabs(a.y), fabs(a.z), fabs(a.w) ); }

double vdot(const DOUBLE4 &a, const DOUBLE4 &b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
}

DOUBLE4 vnormalize(const DOUBLE4 &a)
{
  double l = sqrt(vdot(a, a));
  if (l == 0)
    return 0;
  else
    return a / l;
}


/*
class FLOATN
{
  public:
    float *vm;
    int    vsize;

  FLOATN();
  FLOATN( int vector_size );
  FLOATN( int vector_size, const float *vec );
  FLOATN( const FLOATN &a );
  ~FLOATN();

  void load( int vector_size );
  void load( int vector_size, const float *vec );

  FLOATN& operator=( const FLOATN &a );

  FLOATN operator+( const FLOATN &a ) const;
  FLOATN operator-( const FLOATN &a ) const;
  FLOATN operator*( const FLOATN &a ) const;
  FLOATN operator/( const FLOATN &a ) const;

  FLOATN& operator+=( const FLOATN &a );
  FLOATN& operator-=( const FLOATN &a );
  FLOATN& operator*=( const FLOATN &a );
  FLOATN& operator/=( const FLOATN &a );

  float norm() const;

  FLOATN operator-() const;
};

FLOATN::FLOATN()
{
  memset( this, 0, sizeof(FLOATN) );
}

FLOATN::FLOATN( int vector_size )
{
  memset( this, 0, sizeof(FLOATN) );
  load( vector_size );
}

FLOATN::FLOATN( int vector_size, const float *vec )
{
  memset( this, 0, sizeof(FLOATN) );
  load( vector_size, vec );
}

FLOATN::FLOATN( const FLOATN &a )
{
  memset( this, 0, sizeof(FLOATN) );
  load( a.vsize, a.vm );
}

FLOATN::~FLOATN()
{
  if(vm) free(vm);
}

void FLOATN::load( int vector_size )
{
  if(vm) free(vm);
  vsize = vector_size;
  vm = (float*)malloc( vsize * sizeof(float) );
}

void FLOATN::load( int vector_size, const float *vec )
{
  load( vector_size );
  memcpy( vm, vec, vsize*sizeof(float) );
}

FLOATN& FLOATN::operator=( const FLOATN &a )
{
  load( a.vsize, a.vm );
  return *this;
}

FLOATN FLOATN::operator+( const FLOATN &a ) const
{
  FLOATN res;
    res.load( vsize );

  int i;
  for( i=0; i<vsize; i++ )
    res.vm[i] = vm[i] + a.vm[i];

  return res;
}

FLOATN FLOATN::operator-( const FLOATN &a ) const
{
  FLOATN res;
    res.load( vsize );

  int i;
  for( i=0; i<vsize; i++ )
    res.vm[i] = vm[i] - a.vm[i];

  return res;
}

FLOATN FLOATN::operator*( const FLOATN &a ) const
{
  FLOATN res;
    res.load( vsize );

  int i;
  for( i=0; i<vsize; i++ )
    res.vm[i] = vm[i] * a.vm[i];

  return res;
}

FLOATN FLOATN::operator/( const FLOATN &a ) const
{
  FLOATN res;
    res.load( vsize );

  int i;
  for( i=0; i<vsize; i++ )
    res.vm[i] = vm[i] / a.vm[i];

  return res;
}

FLOATN& FLOATN::operator+=( const FLOATN &a )
{
  int i;
  for( i=0; i<vsize; i++ )
    vm[i] += a.vm[i];

  return *this;
}

FLOATN& FLOATN::operator-=( const FLOATN &a )
{
  int i;
  for( i=0; i<vsize; i++ )
    vm[i] -= a.vm[i];

  return *this;
}

FLOATN& FLOATN::operator*=( const FLOATN &a )
{
  int i;
  for( i=0; i<vsize; i++ )
    vm[i] *= a.vm[i];

  return *this;
}

FLOATN& FLOATN::operator/=( const FLOATN &a )
{
  int i;
  for( i=0; i<vsize; i++ )
    vm[i] /= a.vm[i];

  return *this;
}

float FLOATN::norm() const
{
  return sqrtf( vdot( vm, vm, vsize ) );
}

FLOATN FLOATN::operator-() const
{
  FLOATN res;
    res.load( vsize );

  int i;
  for( i=0; i<vsize; i++ )
    res.vm[i] = -vm[i];

  return res;
}

void YUV2RGB(const FLOAT3 *tYUV, FLOAT3 *tRGB )
{
  float y = tYUV->x;
  float u = tYUV->y;
  float v = tYUV->z;

  float r = y + 1.40200f*u;
  float g = y - 0.71414f*u - 0.34414f*v;
  float b = y              + 1.77200f*v;
  
  tRGB->x = r;
  tRGB->y = g; 
  tRGB->z = b; 

}
void RGB2YUV(FLOAT3 *tYUV, const FLOAT3 *tRGB )
{
  float r = tRGB->x;
  float g = tRGB->y;
  float b = tRGB->z;

  tYUV->x =   0.2990f * r + 0.5870f * g + 0.1140f * b;
  tYUV->y =   0.5000f * r - 0.4187f * g - 0.0813f * b;
  tYUV->z = - 0.1687f * r - 0.3313f * g + 0.5000f * b;
}
*/