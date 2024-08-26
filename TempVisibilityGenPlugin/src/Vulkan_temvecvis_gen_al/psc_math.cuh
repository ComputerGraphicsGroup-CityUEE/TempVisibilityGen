#ifndef PSC_MATH_CUH
#define PSC_MATH_CUH

#include "vec_math.h"

//template <typename T> __device__ T g_clamp( const T &x, const T &a, const T &b){  return x<a ? a : (x<b?x:b);  }
template <typename T> __device__ inline void g_swap( T &a, T &b ){ T c;  c = a;  a = b;  b = c; }

template <typename T> 
__device__ void cg_quick_sort( const T *val, int *idx, int n, bool ascending=true )
{
  if( n<=0 )
    return;

  int low, high, l, h;
  T sv;

  for( l=0; l<n; l++ )
    idx[l] = l;

  int nbuf =   2 * int(log2(double(n)))  +  2;
  int *gs = (int*) malloc( nbuf * sizeof(int) );

  int ns=0;

  gs[ns++]=0;
  gs[ns++]=n-1;
  while( ns )
  {
    high=gs[--ns];
    low=gs[--ns];
    sv = val[idx[high]];

    if(ascending)
    {
      for( l=low, h=low; h<high; h++ )
        if( val[idx[h]]<sv || (val[idx[h]]==sv && idx[h]<idx[high]) )
          g_swap( idx[l++], idx[h] );
    }else
    {
      for( l=low, h=low; h<high; h++ )
        if( val[idx[h]]>sv || (val[idx[h]]==sv && idx[h]<idx[high]) )
          g_swap( idx[l++], idx[h] );
    }
      g_swap( idx[l], idx[h] );
      
    if(ns+4>=nbuf)
    {
      nbuf *= 2;
      int *t = (int*)malloc( nbuf * sizeof(int) );
      memcpy( t, gs, ns*sizeof(int) );
      g_swap( t, gs );
      free(t);
    }
    if( low < l-1 )
    {
      gs[ns++]=low;
      gs[ns++]=l-1;
    }
    if( l+1 < high )
    {
      gs[ns++]=l+1;
      gs[ns++]=high;
    }
  }
  free(gs);
}

template <typename T> 
__device__ void cg_quick_sort( T *val, int n, bool ascending=true )
{
  if( n<=0 )
    return;

  int low, high, l, h;
  float sv;

  //int nbuf =   2 * int(log2(double(n)))  +  2;
  //int *gs = (int*) malloc( nbuf * sizeof(int) );
  int gs[128];

  int ns=0;

  gs[ns++]=0;
  gs[ns++]=n-1;
  while( ns )
  {
    high=gs[--ns];
    low=gs[--ns];
    sv = *((float*)&val[high]);

    if(ascending)
    {
      for( l=low, h=low; h<high; h++ )
        if( *((float*)&val[h])<sv )
          g_swap( val[l++], val[h] );
    }else
    {
      for( l=low, h=low; h<high; h++ )
        if( *((float*)&val[h])>sv )
          g_swap( val[l++], val[h] );
    }
      g_swap( val[l], val[h] );
      
    //if(ns+4>=nbuf)
    //{
    //  nbuf *= 2;
    //  int *t = (int*)malloc( nbuf * sizeof(int) );
    //  memcpy( t, gs, ns*sizeof(int) );
    //  g_swap( t, gs );
    //  free(t);
    //}
    if( low < l-1 )
    {
      if(ns+2>128)
      {
        printf( "cg_quick_sort(), buffer overflowed\n" );
      }
      gs[ns++]=low;
      gs[ns++]=l-1;
    }
    if( l+1 < high )
    {
      if(ns+2>128)
      {
        printf( "cg_quick_sort(), buffer overflowed\n" );
      }
      gs[ns++]=l+1;
      gs[ns++]=high;
    }
  }
  //free(gs);
}


#endif 


