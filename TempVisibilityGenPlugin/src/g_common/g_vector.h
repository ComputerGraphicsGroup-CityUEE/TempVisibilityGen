#ifndef G_VECTOR_H
#define G_VECTOR_H


#define GBYTE unsigned char
class BYTE3;
class BYTE4;
class FLOAT4;


class FLOAT3
{
  public:
    float x;
    float y;
    float z;

  FLOAT3( const BYTE3 &a );
  FLOAT3( float a, float b, float c ):x(a),y(b),z(c){};
  FLOAT3( float a ):x(a),y(a),z(a){};
  FLOAT3():x(0),y(0),z(0){};

  FLOAT3 operator+( const FLOAT3 &a ) const;
  FLOAT3 operator-( const FLOAT3 &a ) const;
  FLOAT3 operator*( const FLOAT3 &a ) const;
  FLOAT3 operator/( const FLOAT3 &a ) const;

  FLOAT3 operator+( const float &f ) const;
  FLOAT3 operator-( const float &f ) const;
  FLOAT3 operator*( const float &f ) const;
  FLOAT3 operator/( const float &f ) const;

  FLOAT3& operator+=( const FLOAT3 &a );
  FLOAT3& operator-=( const FLOAT3 &a );
  FLOAT3& operator*=( const FLOAT3 &a );
  FLOAT3& operator/=( const FLOAT3 &a );

  friend FLOAT3 operator+( float a, const FLOAT3 &b );
  friend FLOAT3 operator-( float a, const FLOAT3 &b );
  friend FLOAT3 operator*( float a, const FLOAT3 &b );
  friend FLOAT3 operator/( float a, const FLOAT3 &b );

  FLOAT3 operator-() const;

  FLOAT3& operator=( float a );

  // row vector affine transform, 
  FLOAT3 rmul( const float *m ) const;
  float norm() const;

  FLOAT3& operator=( const BYTE3 &a );
  bool operator==( const FLOAT3 &b ) const;
  bool operator!=( const FLOAT3 &b ) const;
};

void Angle2Float3( const float theta, const float phi, FLOAT3 *v);
void Float32Angle( float *theta, float *phi, const FLOAT3 *v);

float vdot(const FLOAT3 &a, const FLOAT3 &b);
float vdot(const FLOAT3 *a, const FLOAT3 *b);

FLOAT3 vcross( const FLOAT3 &a, const FLOAT3 &b );

FLOAT3 vnormalize( const FLOAT3 &a );
void vnormalize( FLOAT3 *a );


FLOAT3 vbisect( const FLOAT3 &a, const FLOAT3 &b );
FLOAT3 proj_onto_plane( const FLOAT3 &v, const FLOAT3 &n, const FLOAT3 &v0 );
FLOAT3 intersect( const FLOAT3 &u0, const FLOAT3 &u1, const FLOAT3 &v0, const FLOAT3 &v1 );

FLOAT3 vmax(const FLOAT3 &a, const FLOAT3 &b);
FLOAT3 vmin(const FLOAT3 &a, const FLOAT3 &b);
FLOAT3 sqrtf3( const FLOAT3 &a );
FLOAT3 f3abs( const FLOAT3 &a );
float vangle( const FLOAT3 &a, const FLOAT3 &b );

//void RGB2YUV(       FLOAT3 *tYUV, const FLOAT3 *tRGB );
//void YUV2RGB( const FLOAT3 *tYUV,       FLOAT3 *tRGB );

FLOAT3 MouseOnSphere( int x0, int y0, int width, int height );


void Mul_M33xC3( const float *m, const FLOAT3 *tX, FLOAT3 *tY );
void Mul_R3xM33( const float *m, const FLOAT3 *tX, FLOAT3 *tY );

float vdot ( const float *a, const float *b, int n );
float vnorm( const float *a, int vsize );
float vnormalize( float *a, int n );
void  vsub( const float *a, const float *b, float *c, int vsize );
void  vperturb( const float *a, float *b, int vsize, float perturb_factor );

class FLOAT4
{
  public:
    float x;
    float y;
    float z;
    float w;

  FLOAT4( const FLOAT3 &a, float d ):x(a.x),y(a.y),z(a.z),w(d){};
  FLOAT4( const BYTE4 &a );
  FLOAT4( float a, float b, float c, float d ):x(a),y(b),z(c),w(d){};
  FLOAT4( float a ):x(a),y(a),z(a),w(a){};
  FLOAT4():x(0),y(0),z(0),w(0){};

  FLOAT4 operator+( const FLOAT4 &a ) const;
  FLOAT4 operator-( const FLOAT4 &a ) const;
  FLOAT4 operator*( const FLOAT4 &a ) const;
  FLOAT4 operator/( const FLOAT4 &a ) const;

  FLOAT4 operator+( const float &f ) const;
  FLOAT4 operator-( const float &f ) const;
  FLOAT4 operator*( const float &f ) const;
  FLOAT4 operator/( const float &f ) const;

  FLOAT4& operator+=( const FLOAT4 &a );
  FLOAT4& operator-=( const FLOAT4 &a );
  FLOAT4& operator*=( const FLOAT4 &a );
  FLOAT4& operator/=( const FLOAT4 &a );
  FLOAT4& operator=( const BYTE4 &a );

  friend FLOAT4 operator+( float a, const FLOAT4 &b );
  friend FLOAT4 operator*( float a, const FLOAT4 &b );
  friend FLOAT4 operator-( float a, const FLOAT4 &b );
  friend FLOAT4 operator/( float a, const FLOAT4 &b );

  FLOAT4 operator-() const;

  bool operator==( const FLOAT4 &b ) const;
  bool operator!=( const FLOAT4 &b ) const;
  FLOAT4 operator*( const float *mx ) const;
  friend FLOAT4 operator*( const float *mx, const FLOAT4 &b );
  FLOAT4 rmul( const float *m ) const;
};

FLOAT4 vmax(const FLOAT4 &a, const FLOAT4 &b);
FLOAT4 vmin(const FLOAT4 &a, const FLOAT4 &b);


class BYTE3
{
  public:
    GBYTE x;
    GBYTE y;
    GBYTE z;

  BYTE3& operator=( const GBYTE &a );
  BYTE3& operator=( const FLOAT3 &a );
  bool operator==( const BYTE3 &b ) const;
};

class BYTE4
{
  public:
    GBYTE x;
    GBYTE y;
    GBYTE z;
    GBYTE w;
};




//////////////////////////////////////////////////////////////
// Quaternion Class
//////////////////////////////////////////////////////////////
class GQuat
{
  public:
    float x;
    float y;
    float z;
    float w;

  GQuat():x(0),y(0),z(0),w(1){};
  GQuat(float qx, float qy, float qz, float qw):x(qx),y(qy),z(qz),w(qw){};

  // Notice that :
  // if column vector assumed, GQuat( vFrom, vTo )
  // if row vector assumed,    GQuat( vTo, vFrom )
  GQuat( const FLOAT3 &v0, const FLOAT3 &v1 );

  // Notice that :
  // if column vector assumed, GQuat( axis,  angle )
  // if row vector assumed,    GQuat( axis, -angle )
  GQuat( const FLOAT3 &axis, const float &angle );

  // Quaternion Conjugate
  GQuat operator~() const;

  // Sum of Quaternion
  GQuat operator+( const GQuat& b ) const;

  // Product of Quaternion
  // Caution : Quaternion multiplication DO NOT COMMUTE, therefore
  //           the order of multiplication is IMPORTANT
  GQuat operator*( const GQuat& b ) const;

  // Division of Quaternion
  GQuat operator/( const GQuat& b ) const;

  // Scalar division of Quaternion
  GQuat operator/( const float& b ) const;

  // Scalar multiplication of Quaternion
  GQuat operator*( const float& b ) const;

  // 
  bool operator==( const GQuat& b ) const;

  // Return a normalized version of this Quaternion
  GQuat normalize() const;

  // Filling parameter m with 4x4 rotation matrix
  void matrix( float *m ) const;

  // convert to axis angle representation
  void axis_angle( FLOAT3 &axis, float &angle ) const;

  // Norm of Quaternion
  float norm() const;

  // linear interpolation of quaternions
  static GQuat slerp( const GQuat &a, const GQuat &b, float t );
};

//////////////////////////////////////////////////////////////
// Dint Class
//////////////////////////////////////////////////////////////
class Dnt
{
  public:
    int n;

  Dnt( int bound ):n(bound){};

  int nt( const int &a ){ return (a%n+n)%n; }

};





class DOUBLE3
{
public:
  double x;
  double y;
  double z;

  DOUBLE3(const FLOAT3 &a) :x(a.x), y(a.y), z(a.z) {};
  DOUBLE3(double a, double b, double c) :x(a), y(b), z(c) {};
  DOUBLE3(double a) :x(a), y(a), z(a) {};
  DOUBLE3() :x(0), y(0), z(0) {};

  DOUBLE3 operator+(const DOUBLE3 &a) const;
  DOUBLE3 operator-(const DOUBLE3 &a) const;
  DOUBLE3 operator*(const DOUBLE3 &a) const;
  DOUBLE3 operator/(const DOUBLE3 &a) const;

  DOUBLE3 operator+(const double &f) const;
  DOUBLE3 operator-(const double &f) const;
  DOUBLE3 operator*(const double &f) const;
  DOUBLE3 operator/(const double &f) const;

  DOUBLE3& operator+=(const DOUBLE3 &a);
  DOUBLE3& operator-=(const DOUBLE3 &a);
  DOUBLE3& operator*=(const DOUBLE3 &a);
  DOUBLE3& operator/=(const DOUBLE3 &a);

  friend DOUBLE3 operator+(double a, const DOUBLE3 &b);
  friend DOUBLE3 operator-(double a, const DOUBLE3 &b);
  friend DOUBLE3 operator*(double a, const DOUBLE3 &b);
  friend DOUBLE3 operator/(double a, const DOUBLE3 &b);

  DOUBLE3 operator-() const;
  DOUBLE3 operator+() const;

  DOUBLE3& operator=(double a);
  DOUBLE3& operator=(const FLOAT3 &a);

  // row vector affine transform, 
  DOUBLE3 rmul(const double *m) const;
  double norm() const;

  bool operator==(const DOUBLE3 &b) const;
  bool operator!=(const DOUBLE3 &b) const;
};

double vdot(const DOUBLE3 &a, const DOUBLE3 &b);
DOUBLE3 vnormalize(const DOUBLE3 &a);
DOUBLE3 vcross(const DOUBLE3 &a, const DOUBLE3 &b);
void Angle2Double3(const double theta, const double phi, DOUBLE3 *v);	//convert a solid angle to a vector
void Double32Angle(double *theta, double *phi, const DOUBLE3 *v);

double fmax( const DOUBLE3 &a );


class DOUBLE4
{
public:
  double x;
  double y;
  double z;
  double w;

  DOUBLE4(const DOUBLE4 &a) :x(a.x), y(a.y), z(a.z), w(a.w) {};
  DOUBLE4(double a, double b, double c, double d) :x(a), y(b), z(c), w(d) {};
  DOUBLE4(double a) :x(a), y(a), z(a), w(a) {};
  DOUBLE4() :x(0), y(0), z(0), w(0) {};

  DOUBLE4 operator+(const DOUBLE4 &a) const;
  DOUBLE4 operator-(const DOUBLE4 &a) const;
  DOUBLE4 operator*(const DOUBLE4 &a) const;
  DOUBLE4 operator/(const DOUBLE4 &a) const;

  DOUBLE4 operator+(const double &f) const;
  DOUBLE4 operator-(const double &f) const;
  DOUBLE4 operator*(const double &f) const;
  DOUBLE4 operator/(const double &f) const;

  DOUBLE4& operator+=(const DOUBLE4 &a);
  DOUBLE4& operator-=(const DOUBLE4 &a);
  DOUBLE4& operator*=(const DOUBLE4 &a);
  DOUBLE4& operator/=(const DOUBLE4 &a);

  friend DOUBLE4 operator+(double a, const DOUBLE4 &b);
  friend DOUBLE4 operator-(double a, const DOUBLE4 &b);
  friend DOUBLE4 operator*(double a, const DOUBLE4 &b);
  friend DOUBLE4 operator/(double a, const DOUBLE4 &b);
  friend DOUBLE4 operator*(double* m, const DOUBLE4 &b);

  DOUBLE4 operator-() const;
  DOUBLE4 operator+() const;

  DOUBLE4& operator=(double a);

  double norm() const;

  bool operator==(const DOUBLE4 &b) const;
  bool operator!=(const DOUBLE4 &b) const;
};

DOUBLE4 log( const DOUBLE4 &a );
DOUBLE4 tanh( const DOUBLE4 &a );
DOUBLE4 fabs( const DOUBLE4 &a );
double vdot(const DOUBLE4 &a, const DOUBLE4 &b);
DOUBLE4 vnormalize(const DOUBLE4 &a);















#endif