//#include <stdlib.h> 
//#include <time.h>

#include "gen_vectorized_visibility.h"
#include "bvh_cuda.h"
#include "cal_visibility_cuda.h"
#include "g_common.h"


void save_PointIArr(char *spath, PointI* pibuf, int npi)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<npi; i++)
  {
		fprintf(f0, "%d\t%f\t",pibuf[i].ei,pibuf[i].es);
    if((i+1)%2 ==0)
      fprintf(f0, "\n");
  }

	fclose(f0);
}

void save_intArr(char *spath, int *ibuf, int nbuf, int npl)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nbuf; i++)
	{
		fprintf(f0, "%d\t",ibuf[i]);
		if((i+1)%npl==0)
			fprintf(f0, "\n");
	}

	fclose(f0);
}

#include "g_pfm.h"
void save_PsArr(char *spath, cups *ibuf, int nbuf)
{
  GPf1 idx0, idx1;
  GPfm emap0, emap1, omap, vnmap;
  int w = 512;
  int h = nbuf/w +1;
  idx0.load(w, h);
  idx1.load(w, h);
  emap0.load(w, h);
  emap1.load(w, h);
  omap.load(w, h);
  vnmap.load(w, h);
  printf("nbuf, w, h %d %d %d\n", nbuf, w, h);

  for(int i=0; i<w*h;i++)
  {
    if(i<nbuf)
    {
      idx0.fm[i] = ibuf[i].eidx;
      idx1.fm[i] = ibuf[i].tri_idx;
      emap0.fm[i] = FLOAT3(ibuf[i].ev0.x,ibuf[i].ev0.y,ibuf[i].ev0.z);
      emap1.fm[i] =  FLOAT3(ibuf[i].ev1.x,ibuf[i].ev1.y,ibuf[i].ev1.z);
      omap.fm[i] = FLOAT3(ibuf[i].origin.x,ibuf[i].origin.y,ibuf[i].origin.z);
      vnmap.fm[i] = FLOAT3(ibuf[i].vvn.x,ibuf[i].vvn.y,ibuf[i].vvn.z);
     }
  }

  idx0.save("idx0_dll.pfm");
  idx1.save("idx1_dll.pfm");
  emap0.save("emap0_dll.pfm");
  emap1.save("emap1_dll.pfm");
  omap.save("omap_dll.pfm");
  vnmap.save("vnmap_dll.pfm");


}

void save_intPf1(char *spath, int *ibuf, int nibuf)
{
  int w = 512;
  int h = nibuf/w +1;
  GPf1 dst;
  dst.load(w, h);

  for(int i=0; i<w*h; i++)
  {
    if(i<nibuf)
      dst.fm[i] = ibuf[i];
  }
  dst.save(spath);
}



void save_cubvh(cubvh* bvhbuf, int nb)
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
void save_cubvhArr(std::vector<cubvh> &bvhbuf, int nb)
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



VisGen::VisGen()
{
  memset(this,0,sizeof(VisGen));
  max_front_back_edge_buf = 256;
  total_line_sample = 0;
}

VisGen::~VisGen()
{
  free_cuda_memory();
}

void VisGen::init_cuda_memory(
  int n_vertex,
  int *face_vidx, int *face_eidx, int n_face,
  GEdge* s_edge, int n_edge
){
  //printf( ">>> vmap %f MB\n", float(n_vertex*sizeof(float3))/1024/1024 );
  //printf( ">>> nmap %f MB\n", float(n_face*sizeof(float3))/1024/1024 );
  //printf( ">>> emap %f MB\n", float(n_edge*sizeof(cuedge))/1024/1024 );
  //printf( ">>> vimap %f MB\n", float(n_face*3*sizeof(int))/1024/1024 );
  //printf( ">>> eimap %f MB\n", float(n_face*3*sizeof(int))/1024/1024 );
  //printf( ">>> vnmap %f MB\n", float(n_vertex*sizeof(float3))/1024/1024 );
  cudaMalloc( (void**)&vis_param.vmap, n_vertex*sizeof(float3) );
  cudaMalloc( (void**)&vis_param.nmap, n_face*sizeof(float3) );
  cudaMalloc( (void**)&vis_param.emap, n_edge*sizeof(cuedge) );
  cudaMalloc( (void**)&vis_param.vimap, n_face*3*sizeof(int) );
  cudaMalloc( (void**)&vis_param.eimap, n_face*3*sizeof(int) );
  cudaMalloc( (void**)&vis_param.vnmap, n_vertex*sizeof(float3) );


  cudaMemcpy(vis_param.emap, s_edge, n_edge*sizeof(GEdge), cudaMemcpyHostToDevice); // sizeof(Gedge)
  cudaMemcpy(vis_param.vimap, face_vidx, n_face*3*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vis_param.eimap, face_eidx, n_face*3*sizeof(int), cudaMemcpyHostToDevice);

  int max_vertex_perpass;
  max_vertex_perpass = floor( floor( double(max_front_back_edge_buf) *1024 *1024 /n_edge /4 )/4 )*4;
  n_vertex_perpass = g_min( n_vertex, max_vertex_perpass );
  //printf(">>> n_vertex_perpass %i\n", n_vertex_perpass);

  //printf( ">>> front_back_edge_buf %f MB\n", float(n_vertex_perpass*n_edge*sizeof(int))/1024/1024 );
  //printf( ">>> front_back_edge_offset %f MB\n", float((((n_edge+(32-1))/32) * n_vertex_perpass  + 1 )*sizeof(int))/1024/1024 );
  //printf( ">>> front_back_edge_vinfo %f MB\n", float((n_vertex+1)*sizeof(int))/1024/1024 );
  //printf( ">>> vis_tightpacked_vinfo %f MB\n", float((n_vertex+1)*sizeof(int))/1024/1024 );
  cudaMalloc( (void**)&front_back_edge_buf, n_vertex_perpass*n_edge*sizeof(int) );
  cudaMalloc( (void**)&front_back_edge_offset, (((n_edge+(32-1))/32) * n_vertex_perpass  + 1 )*sizeof(int) );
  cudaMemset( front_back_edge_offset, 0, sizeof(int));
  cudaMalloc((void**)&front_back_edge_vinfo, (n_vertex+1)*sizeof(int));
  cudaMalloc((void**)&vis_tightpacked_vinfo, (n_vertex+1)*sizeof(int));

}

void VisGen::free_cuda_memory()
{
  cudaFree(vis_param.cbvh);
  cudaFree(vis_param.vmap);
  cudaFree(vis_param.nmap);
  cudaFree(vis_param.emap);
  cudaFree(vis_param.vimap);
  cudaFree(vis_param.eimap);

  cudaFree(vis_param.vnmap);

  cudaFree(vis_record);

  cudaFree(front_back_edge_buf);
  cudaFree(front_back_edge_offset);
  cudaFree(front_back_edge_vinfo);
  cudaFree(vis_tightpacked_vinfo);

  cudaFree(front_back_edge_tightpacked);
  cudaFree(ps_buf);
  cudaFree(vis_buf);
  cudaFree(vis_offset);
  cudaFree(vis_tightpacked);


  if(my_dat_buf)
    free(my_dat_buf);
  if(my_dat_vi)
    free(my_dat_vi);

  memset(this,0,sizeof(VisGen));
}

void VisGen::update_cuda_memory(
  MyVertex *vbuf0, const std::vector< std::vector<uint32_t> > &vidx_distinct_all,
  FLOAT3 *fnbuf, int n_face, int n_edge
){  
  int n_vertex = int(vidx_distinct_all.size());

  CUParam &param = vis_param;

  FLOAT3 *vbuf, *vnbuf;
  vbuf  = (FLOAT3*) malloc( n_vertex*sizeof(FLOAT3) );
  vnbuf = (FLOAT3*) malloc( n_vertex*sizeof(FLOAT3) );
  int i;
  for( i=0; i<n_vertex; i++ )
  {
    vbuf[i]  = *((FLOAT3*)&vbuf0[ vidx_distinct_all[i][0] ].pos);
    vnbuf[i] = *((FLOAT3*)&vbuf0[ vidx_distinct_all[i][0] ].normal);
  }

  int max_vertex_perpass;
  max_vertex_perpass = floor( floor( double(max_front_back_edge_buf) *1024 *1024 /n_edge /4 )/4 )*4;
  n_vertex_perpass = g_min( n_vertex, max_vertex_perpass );

  //printf("update_cuda_memory n_vertex_perpass %i\n", n_vertex_perpass);


  cudaMemcpy(param.nmap, fnbuf, n_face*sizeof(FLOAT3), cudaMemcpyHostToDevice);
  cudaMemcpy(param.vmap, vbuf, n_vertex*sizeof(FLOAT3), cudaMemcpyHostToDevice);
  cudaMemcpy(param.vnmap, vnbuf, n_vertex*sizeof(FLOAT3), cudaMemcpyHostToDevice);

  free(vbuf);
  free(vnbuf);
}


void VisGen::update_cuda_memoryX(
  MyVertex *vbuf0, const std::vector<int> &vidx_distinct_all,
  FLOAT3 *fnbuf, int n_face, int n_edge
){  
  int n_vertex = int(vidx_distinct_all.size());

  CUParam &param = vis_param;

  FLOAT3 *vbuf, *vnbuf;
  vbuf  = (FLOAT3*) malloc( n_vertex*sizeof(FLOAT3) );
  vnbuf = (FLOAT3*) malloc( n_vertex*sizeof(FLOAT3) );
  int i;
  for( i=0; i<n_vertex; i++ )
  {
    vbuf[i]  = *((FLOAT3*)&vbuf0[ vidx_distinct_all[i] ].pos);
    vnbuf[i] = *((FLOAT3*)&vbuf0[ vidx_distinct_all[i] ].normal);
  }

  int max_vertex_perpass;
  max_vertex_perpass = floor( floor( double(max_front_back_edge_buf) *1024 *1024 /n_edge /4 )/4 )*4;
  n_vertex_perpass = g_min( n_vertex, max_vertex_perpass );

  //printf("update_cuda_memory n_vertex_perpass %i\n", n_vertex_perpass);

  cudaMemcpy(param.nmap, fnbuf, n_face*sizeof(FLOAT3), cudaMemcpyHostToDevice);
  cudaMemcpy(param.vmap, vbuf, n_vertex*sizeof(FLOAT3), cudaMemcpyHostToDevice);
  cudaMemcpy(param.vnmap, vnbuf, n_vertex*sizeof(FLOAT3), cudaMemcpyHostToDevice);

  free(vbuf);
  free(vnbuf);
}




void VisGen::update_bvh( cubvh *dat, int node_count )
{
  CUParam &param = vis_param;
  //cubvh *dat;
  //bvh2cuda(bvh0, dat);

  vis_bvh0_node_count = node_count;
  if( _vis_bvh0_node_count < vis_bvh0_node_count )
  {
    _vis_bvh0_node_count = int(1.5*vis_bvh0_node_count);
    cudaFree(param.cbvh);

    //printf( ">>> cbvh %f MB\n", float(_vis_bvh0_node_count*sizeof(cubvh))/1024/1024 );
    cudaMalloc( (void**)&param.cbvh, _vis_bvh0_node_count*sizeof(cubvh) );

    //printf( "VisGen::gen_vis_rendering3(), memory allocated %d\n", _vis_bvh0_node_count);
  }


  cudaMemcpy(param.cbvh, dat, node_count*sizeof(cubvh), cudaMemcpyHostToDevice);
  //save_cubvh(dat, node_count);
  //free(dat);
}


void VisGen::update_bvh(std::vector<cubvh> &dat, int node_count)
{
  CUParam &param = vis_param;
  //cubvh *dat;
  //bvh2cuda(bvh0, dat);

  vis_bvh0_node_count = node_count;
  if( _vis_bvh0_node_count < vis_bvh0_node_count )
  {
    _vis_bvh0_node_count = int(1.5*vis_bvh0_node_count);
    cudaFree(param.cbvh);

    printf( ">>> cbvh %f MB\n", float(_vis_bvh0_node_count*sizeof(cubvh))/1024/1024 );
    cudaMalloc( (void**)&param.cbvh, _vis_bvh0_node_count*sizeof(cubvh) );

    printf( "VisGen::gen_vis_rendering3(), memory allocated %d\n", _vis_bvh0_node_count);
  }


  cudaMemcpy(param.cbvh, dat.data(), node_count*sizeof(cubvh), cudaMemcpyHostToDevice);
  save_cubvhArr(dat, node_count);
}



void VisGen::gen_vis_rendering3(
  int n_edge, int vertex0, int n_vertex
){
  CUParam &param = vis_param;

  //my_check_cuda();
  //printf("front_back_edge_tightpacked_size %d\n", front_back_edge_tightpacked_size);
  //printf("_front_back_edge_tightpacked_size %d\n\n", _front_back_edge_tightpacked_size);
  //printf("n_edge %d\n\n", n_edge);
  //printf("vertex0 %d\n\n", vertex0);
  //printf("n_vertex %d\n\n", n_vertex);


  cu_initialize_ps(
    front_back_edge_buf, 
    front_back_edge_offset, 
    ps_buf, 
    front_back_edge_tightpacked, 
    front_back_edge_tightpacked_size, 
   _front_back_edge_tightpacked_size, 
    front_back_edge_vinfo, 
    param, n_edge, 
    
    vertex0,
    n_vertex
  
  );

  //printf("n_edge %d\n\n", n_edge);
  //printf("vertex0 %d\n\n", vertex0);
  //printf("n_vertex %d\n\n", n_vertex);

  //my_check_cuda();
  {
  //cups* tmppsbuf;
  // tmppsbuf =  (cups*) malloc((_front_back_edge_tightpacked_size+(n_vertex*32))*sizeof(cups));
  //cudaMemcpy( tmppsbuf, ps_buf,(_front_back_edge_tightpacked_size+(n_vertex*32))*sizeof(cups), cudaMemcpyDeviceToHost );
  //save_PsArr("ps_buf_dll.txt", tmppsbuf, _front_back_edge_tightpacked_size+(n_vertex*32));
  //free(tmppsbuf);
  //int *tmpfbbuf;
  //printf("n_vertex_perpass, n_edge %d %d\n",n_vertex_perpass, n_edge);
  //tmpfbbuf = (int*)malloc(n_vertex_perpass*n_edge*sizeof(int));
  //cudaMemcpy( tmpfbbuf, front_back_edge_buf, n_vertex_perpass*n_edge*sizeof(int), cudaMemcpyDeviceToHost );
  //save_intPf1("front_back_edge_buf_dll.pfm", tmpfbbuf, n_vertex_perpass*n_edge);
  //free(tmpfbbuf);
  //int *tmpfbvinfo;
  //tmpfbvinfo = (int*)malloc((n_vertex+1)*sizeof(int));
  //cudaMemcpy( tmpfbvinfo, front_back_edge_vinfo, (n_vertex+1)*sizeof(int), cudaMemcpyDeviceToHost );
  //save_intPf1("front_back_edge_vinfo_dll.pfm", tmpfbvinfo, (n_vertex+1));
  //free(tmpfbvinfo);

  //int *tmpfboff;
  //tmpfboff = (int*)malloc((((n_edge+(32-1))/32) * n_vertex_perpass  + 1 )*sizeof(int));
  //cudaMemcpy( tmpfboff, front_back_edge_offset, (((n_edge+(32-1))/32) * n_vertex_perpass  + 1 )*sizeof(int), cudaMemcpyDeviceToHost );
  //save_intPf1("front_back_edge_offset_dll.pfm", tmpfboff, (((n_edge+(32-1))/32) * n_vertex_perpass  + 1 ));
  //free(tmpfboff);

  //int *tmpfbtight;
  //tmpfbtight = (int*)malloc(2*_front_back_edge_tightpacked_size*sizeof(int));
  //cudaMemcpy( tmpfbtight, front_back_edge_tightpacked, 2*_front_back_edge_tightpacked_size * sizeof(int), cudaMemcpyDeviceToHost );
  //save_intPf1("front_back_edge_tightpacked_dll.pfm", tmpfbtight, 2*_front_back_edge_tightpacked_size);
  //free(tmpfbtight);

  //printf("front_back_edge_tightpacked_size %d\n", front_back_edge_tightpacked_size);
  //printf("_front_back_edge_tightpacked_size %d\n", _front_back_edge_tightpacked_size);

  //save_intArrX("front_back_edge_vinfo.txt", front_back_edge_vinfo,n_vertex+1, 4);
  }
  nps = front_back_edge_tightpacked_size + n_vertex*32;

  total_line_sample += nps;
  //printf("nps, total_line_sample %d %d\n", nps, total_line_sample);

  if( _nps<nps )
  {
    _nps = int(nps*1.5);
    if(vis_buf)
      cudaFree(vis_buf);
    if(vis_offset)
      cudaFree(vis_offset);

    //printf( ">>> vis_buf %f MB\n", float(_nps)*RECORD_EDGE_SIZE*sizeof(PointI)/1024/1024 );
    //printf( ">>> vis_offset %f MB\n", float((_nps+1)*sizeof(int))/1024/1024 );
    cudaMalloc( (void**)&vis_buf, _nps*RECORD_EDGE_SIZE*sizeof(PointI) );
    cudaMalloc((void**)&vis_offset, (_nps+1)*sizeof(int));
    cudaMemset(vis_offset, 0, sizeof(int));

    //printf( "VisGen::gen_vis_rendering3(), memory allocated %s\n", "vis_buf" );
    //printf( "VisGen::gen_vis_rendering3(), memory allocated %s\n", "vis_offset" );
    //my_check_cuda();
  }

  //printf( "VisGen::gen_vis_rendering3(), cu_ps_bvh_flattened_cone_speed_g0()...\n" );
  //{
    //printf("vis_tightpacked_size %d\n", vis_tightpacked_size);
    //printf("_vis_tightpacked_size %d\n", _vis_tightpacked_size);
    //printf("nps %d\n\n", nps);

  //}

  cu_ps_bvh_flattened_cone_speed_g0(
    vis_tightpacked, vis_tightpacked_size, _vis_tightpacked_size, vis_tightpacked_vinfo, 
    vis_buf, vis_offset,
    param, ps_buf, front_back_edge_vinfo, 
    nps, 
    //g_min( n_vertex, n_vertex_perpass ) 
    n_vertex
  );
  //my_check_cuda();
  //printf( "done\n" );

  //{
    //printf("vis_tightpacked_size %d\n", vis_tightpacked_size);
    //printf("_vis_tightpacked_size %d\n", _vis_tightpacked_size);
    //printf("nps %d\n", nps);
  //}

  my_dat_size = vis_tightpacked_size;
  //printf("my_dat_size %d\n", my_dat_size);

  if( _my_dat_size < my_dat_size )
  {
    //printf( "VisGen::gen_vis_rendering3(), my_dat allocate... " );
    if(my_dat_vi)
      free(my_dat_vi);
    if(my_dat_buf)
      free(my_dat_buf);
    _my_dat_size = int(1.5*my_dat_size);
    my_dat_vi = (int*) malloc( (n_vertex+1)*sizeof(int) );
    my_dat_buf = (PointI*) malloc( _my_dat_size*sizeof(PointI) );
    //printf( "done\n" );
  }

  //printf( "VisGen::gen_vis_rendering3(), copy my_dat_buf... " );
  cudaMemcpy( my_dat_buf, vis_tightpacked, vis_tightpacked_size*sizeof(PointI), cudaMemcpyDeviceToHost );
  my_check_cuda();
  //printf( "done\n" );

  //printf( "VisGen::gen_vis_rendering3(), copy my_dat_vi... " );
  cudaMemcpy( my_dat_vi, vis_tightpacked_vinfo, (n_vertex+1)*sizeof(int), cudaMemcpyDeviceToHost );
  my_check_cuda();
  //printf( "done\n" );


  //save_PointIArr("DLL_my_dat_buf.txt", my_dat_buf, vis_tightpacked_size);
  //save_intArr("DLL_my_dat_vi.txt", my_dat_vi,n_vertex+1, 4);
}


void VisGen::testing( int n_edge, int n_vertex, int n_update )
{
  CUParam &param = vis_param;


  cu_initialize_ps(
    front_back_edge_buf, 
    front_back_edge_offset, 
    ps_buf, 
    front_back_edge_tightpacked, 
    front_back_edge_tightpacked_size, 
   _front_back_edge_tightpacked_size, 
    front_back_edge_vinfo, 
    param, n_edge, n_vertex, n_update);

  nps = front_back_edge_tightpacked_size + n_vertex*32;
  //printf("ps %d\n", nps);

  if( _nps<nps )
  {
    _nps = int(nps*1.5);
    if(vis_buf)
      cudaFree(vis_buf);
    if(vis_offset)
      cudaFree(vis_offset);
    printf( ">>> vis_buf %f MB\n", float(_nps*RECORD_EDGE_SIZE*sizeof(PointI))/1024/1024 );
    printf( ">>> vis_offset %f MB\n", float((_nps+1)*sizeof(int))/1024/1024 );
    cudaMalloc( (void**)&vis_buf, _nps*RECORD_EDGE_SIZE*sizeof(PointI) );
    cudaMalloc((void**)&vis_offset, (_nps+1)*sizeof(int));
    cudaMemset(vis_offset, 0, sizeof(int));


    printf( "VisGen::gen_vis_rendering3(), memory allocated %s\n", "vis_buf" );
    printf( "VisGen::gen_vis_rendering3(), memory allocated %s\n", "vis_offset" );
  }

  cu_ps_bvh_flattened_cone_speed_testing(
    vis_buf, vis_offset,
    param, ps_buf, front_back_edge_vinfo, nps, n_vertex, n_update
  );

}





