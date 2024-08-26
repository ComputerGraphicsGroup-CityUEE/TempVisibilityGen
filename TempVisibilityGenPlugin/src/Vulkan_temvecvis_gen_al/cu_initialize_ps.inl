__device__ bool isfrontbackedge(const float3 &pp, const float3 &n0, const float3 &n1)
{
  if(dot(pp, n0) > 0 != dot(pp, n1) > 0)
    return true;
  else
    return false;
}

__device__ int determine_ps(CUParam &param, int vidx, int eidx )
{
  float3 my_rand = make_float3(-0.520755f, 0.107965f, -0.846852f);
  float3 origin = param.vmap[vidx];
  float3 vvn = param.vnmap[vidx];
  origin += vvn*.002 + my_rand*.001;


  float3 v0, v1, n0, n1, pp;
  v0 = param.vmap[param.emap[eidx].v0]; 
  v1 = param.vmap[param.emap[eidx].v1];
  n0 = param.nmap[param.emap[eidx].t0];

  if( param.emap[eidx].t1==-1 )
    n1 = -n0;
  else
    n1 = param.nmap[param.emap[eidx].t1];

  pp = origin - v0;
  if(
    (  param.emap[eidx].t1==-1 || isfrontbackedge(pp, n0, n1)  ) 
        && ( dot(v0-origin, vvn)>0 || dot(v1-origin, vvn)>0 )
  ){   
    return eidx+1;
  }else
    return 0;

}

__device__ void cal_prefix_sum( int *prefix_sum1, int *prefix_sum0, int i, int id )
{
  if(i-1<id)
    prefix_sum1[id] = prefix_sum0[id-i] + prefix_sum0[id];
  else
    prefix_sum1[id] =                     prefix_sum0[id];
}


__global__ void call_initi_ps( CUParam param, int *front_back_edge_buf, int *front_back_edge_offset, int n_edge, int vertex0 )
{  
  // blockDim.x 32;  [0,  32) <- threadIdx.x
  // gridDim.x 166;  [0, 166) <- blockIdx.x
  // gridDim.y 1941; [0,1941) <- blockIdx.y

 __shared__ int res[32], prefix_sum0[32], prefix_sum1[32];

  int id = threadIdx.x;
  int e0 = blockIdx.x*blockDim.x;
  int ei = e0 + threadIdx.x;

  int vi = blockIdx.y + vertex0;
  int vi_per_pass = blockIdx.y;


  res[id] = 0;
  prefix_sum0[id] = 0;
  if( ei < n_edge )
  {
    res[id] = determine_ps(param, vi, ei);
    if(res[id])
      prefix_sum0[id] = 1;
  }
  __syncthreads();

  cal_prefix_sum( prefix_sum1, prefix_sum0, 1, id );
  __syncthreads();
  cal_prefix_sum( prefix_sum0, prefix_sum1, 2, id );
  __syncthreads();
  cal_prefix_sum( prefix_sum1, prefix_sum0, 4, id );
  __syncthreads();
  cal_prefix_sum( prefix_sum0, prefix_sum1, 8, id );
  __syncthreads();
  cal_prefix_sum( prefix_sum1, prefix_sum0, 16, id );
  __syncthreads();

  if(res[id])
    front_back_edge_buf[ vi_per_pass*n_edge + e0 + prefix_sum1[id]-1] = res[id]-1;
  if(id==0)
    front_back_edge_offset[ vi_per_pass*gridDim.x + blockIdx.x ] = prefix_sum1[32-1];
}


__global__ void tight_pack2( int *front_back_edge_tightpacked, int *front_back_edge_buf, int *front_back_edge_offset, int n_edge, int E1, int vertex0 )
{  
  // gridDim.x   22;  [0,  22) <- blockIdx.x
  // gridDim.y  328;  [0, 328) <- blockIdx.y
  // blockDim.x 256;  [0, 256) <- threadIdx.x

  int E1id = blockIdx.x * blockDim.x + threadIdx.x;

  if( E1id<E1 )
  {
    int vi_per_pass = blockIdx.y;
    int vi = blockIdx.y + vertex0;

    int id  = vi_per_pass * E1 + E1id;
    int idx = vi_per_pass * n_edge  + E1id*32;
    int i0 = front_back_edge_offset[id];
    int i1 = front_back_edge_offset[id+1];
    for( int i=0; i<i1-i0; i++ )
    {
      front_back_edge_tightpacked[2*(i0+i)+0] = vi;
      front_back_edge_tightpacked[2*(i0+i)+1] = front_back_edge_buf[idx+i];
    }
  }



  //int blockId = gridDim.x * gridDim.y*blockIdx.z + gridDim.x*blockIdx.y + blockIdx.x;
  //int blockSize = blockDim.x;
  //int blockId = gridDim.x*blockIdx.y + blockIdx.x;
  //int threadId = blockId * blockDim.x + threadIdx.x;
  //int E1id = threadId; //blockDim.x*blockIdx.y +threadIdx.x;

  //if( E1id<E1 )
  //{
  //  int id  = blockIdx.x * E1 + E1id;
  //  int idx = blockIdx.x * n_edge  + E1id*32;
  //  int i0 = front_back_edge_offset[id];
  //  int i1 = front_back_edge_offset[id+1];
  //  int vi = blockIdx.x;
  //  for( int i=0; i<i1-i0; i++ )
  //  {
  //    front_back_edge_tightpacked[2*(i0+i)+0] = vi;
  //    front_back_edge_tightpacked[2*(i0+i)+1] = front_back_edge_buf[idx+i];
  //  }
  //}





}


__global__ void update_index( int *front_back_edge_vinfo, int *front_back_edge_offset, int n_block, int n_vertex)
{
  // blockDim.x 256;  [0,  256) <- threadIdx.x
  // gridDim.x (1941+1+(256-1))/256;  [0, (1941+1+(256-1))/256) <- blockIdx.x
  int id = blockIdx.x*blockDim.x + threadIdx.x;

  if(id>=n_vertex)
    return;

  front_back_edge_vinfo[id] = front_back_edge_offset[id*n_block];
}


__global__ void calssss( cups *ps_buf, CUParam param, int *front_back_edge_vinfo, int *front_back_edge_tightpacked, int front_back_edge_tightpacked_size, int vertex0)
{
  // blockDim.x 256;  [0,  256) <- threadIdx.x
  // gridDim.x (front_back_edge_tightpacked_size+(256-1))/256;  [0, (buf_size+(256-1))/256) <- blockIdx.x

  int i1 = blockIdx.x*blockDim.x + threadIdx.x;

  if(i1>=front_back_edge_tightpacked_size)
    return;


  int vi = front_back_edge_tightpacked[2*i1+0];
  int vi_per_pass = vi-vertex0;
  int ei = front_back_edge_tightpacked[2*i1+1];
  int i0 = front_back_edge_vinfo[vi_per_pass];
  int ii = i1-i0;

  float3 my_rand = make_float3(-0.520755f, 0.107965f, -0.846852f);
  cups ps;
    ps.vvn = param.vnmap[vi];
    ps.origin = param.vmap[vi]+ ps.vvn*.002 + my_rand*.001;
    ps.eidx = ei;
    float3 vn = param.nmap[param.emap[ps.eidx].t0];
    if (dot(vn, ps.origin - param.vmap[param.emap[ps.eidx].v0]) > 0)
    {
      ps.tri_idx = param.emap[ps.eidx].t0;
      ps.ev0 = param.vmap[ param.emap[ps.eidx].v0 ];
      ps.ev1 = param.vmap[ param.emap[ps.eidx].v1 ];
    }else
    {
      ps.tri_idx = param.emap[ps.eidx].t1;
      ps.ev0 = param.vmap[ param.emap[ps.eidx].v1 ];
      ps.ev1 = param.vmap[ param.emap[ps.eidx].v0 ];
    }

  ps_buf[  i0+vi_per_pass*32 + ii  ] = ps;
}


__global__ void caltttt( cups *ps_buf, CUParam param, int *front_back_edge_vinfo, int vertex0 )
{

  // blockDim.x 32;  [0,  32) <- threadIdx.x
  // gridDim.x n_vertex;  [0, n_vertex) <- blockIdx.x


  int vi_per_pass = blockIdx.x;
  int vi = vertex0 + vi_per_pass;

  int i0 = front_back_edge_vinfo[vi_per_pass+1] + 32*vi_per_pass;

  int nstep = blockDim.x;
  int ii = threadIdx.x;
  float3 my_rand = make_float3(-0.520755f, 0.107965f, -0.846852f);
  cups ps;
    ps.vvn = param.vnmap[vi];
    ps.origin = param.vmap[vi]+ ps.vvn*.002 + my_rand*.001;
    float3 ax, ay, az;
    get_coordinate_system(ps, ax, ay, az);

    ps.eidx = -1;
    ps.tri_idx = -1;

  float theta0 = float( ii+0       )/nstep*2*G_PI;
  float theta1 = float((ii+1)%nstep)/nstep*2*G_PI;
  ps.ev1 = ps.origin + (cosf(theta0)*ax + sinf(theta0)*az)*2;
  ps.ev0 = ps.origin + (cosf(theta1)*ax + sinf(theta1)*az)*2;

  ps_buf[  i0+ii   ] = ps;
}


__host__ void cu_initialize_ps(
  int *front_back_edge_buf, 
  int *front_back_edge_offset, 
  cups *&ps_buf,
  int *&front_back_edge_tightpacked, 
  int &front_back_edge_tightpacked_size, 
  int &_front_back_edge_tightpacked_size,
  int *front_back_edge_vinfo,
  const CUParam &param, int n_edge, int vertex0, int n_vertex )
{
  //printf( "cu_initialize_ps(): vertex0 %i\n", vertex0 );
  //printf( "cu_initialize_ps(): n_vertex %i\n", n_vertex );
  //printf( "\n" );

  //printf( "cu_initialize_ps(): call_initi_ps... " );
  call_initi_ps<<<  dim3((n_edge+(32-1))/32,
    n_vertex
  ), 32  >>>(param, front_back_edge_buf, front_back_edge_offset+1, n_edge, 
    vertex0
    );
  //printf( "done\n" );
  //my_check_cuda();


  //printf( "cu_initialize_ps(): thrust::inclusive_scan... " );
  //thrust::inclusive_scan(thrust::device, 
  //  front_back_edge_offset+1 + (((n_edge+(32-1))/32) * vertex0-1), 
  //  front_back_edge_offset+1 + ((n_edge+(32-1))/32) * n_vertex, 
  //  front_back_edge_offset+1 + (((n_edge+(32-1))/32) * vertex0-1)
  //);
  thrust::inclusive_scan(thrust::device, 
    front_back_edge_offset,
    front_back_edge_offset+1 + ((n_edge+(32-1))/32) * n_vertex, 
    front_back_edge_offset
  );
  //printf( "done\n" );
  //my_check_cuda();

  cudaMemcpy( &front_back_edge_tightpacked_size, front_back_edge_offset+((n_edge+(32-1))/32)*n_vertex, sizeof(int), cudaMemcpyDeviceToHost );
  //printf( "front_back_edge_tightpacked_size %i\n", front_back_edge_tightpacked_size );
  //my_check_cuda();

  if(_front_back_edge_tightpacked_size<front_back_edge_tightpacked_size)
  {
    _front_back_edge_tightpacked_size = int(1.5*front_back_edge_tightpacked_size);
    if(front_back_edge_tightpacked)
      cudaFree(front_back_edge_tightpacked);
    if(ps_buf)
      cudaFree(ps_buf);
    cudaMalloc( (void**)&front_back_edge_tightpacked, 2*_front_back_edge_tightpacked_size * sizeof(int) );
    cudaMalloc( (void**)&ps_buf, (_front_back_edge_tightpacked_size+(n_vertex*32)) * sizeof(cups) );
    //printf( "cal_initialize_ps(), memory allocated\n" );
    //my_check_cuda();
  }

  //printf( "tight_pack2()...  " );
  int E1 = (n_edge+(32-1))/32;
  int BB = 256;
  tight_pack2<<<  dim3((E1+(BB-1))/BB, n_vertex ), BB  >>>
    ( front_back_edge_tightpacked, front_back_edge_buf, front_back_edge_offset, n_edge,E1, vertex0 );
  //printf( "done\n" );
  //my_check_cuda();

  //printf( "update_index()...  " );
  update_index<<< (n_vertex+1+(256-1))/256, 256 >>>(front_back_edge_vinfo, front_back_edge_offset, (n_edge+(32-1))/32, n_vertex+1);
  //printf( "done\n" );
  //my_check_cuda();

  //printf( "calssss()...  " );
  calssss<<<  (front_back_edge_tightpacked_size +(256-1))/256   , 256  >>>(ps_buf, param, front_back_edge_vinfo, front_back_edge_tightpacked, front_back_edge_tightpacked_size, vertex0);
  //printf( "done\n" );
  //my_check_cuda();

  //printf( "caltttt()...  " );
  caltttt<<<  n_vertex, 32  >>>(ps_buf, param, front_back_edge_vinfo, vertex0 );
  //printf( "done\n" );
  //my_check_cuda();


  //cudaError_t cudaStatus;
  //cudaStatus = cudaGetLastError();
  //if( cudaStatus != cudaSuccess )
  //  printf( "CUDA error: %s\n", cudaGetErrorString(cudaStatus));
  //cudaStatus = cudaDeviceSynchronize();
  //if( cudaStatus != cudaSuccess )
  //  printf( "cudaDeviceSynchronize error: %s\n", cudaGetErrorString(cudaStatus));
}





