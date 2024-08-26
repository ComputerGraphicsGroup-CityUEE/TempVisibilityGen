// dllmain.cpp : Defines the entry point for the DLL application.
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <Windows.h>
#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"

#include "shaders/my_raytrace_datatypes.h"
#include "bvh_test.h"
#include "cal_visibility_cuda.h"
#include "g_common.h"

#include "Qtype.h"
#include "genvis_assis_func.h"

extern "C" __declspec(dllexport) void __stdcall
gen_tempVecVis_ai(
  const Vector3 *_ls_v, 
  const Vector3 *_ls_vn, // for vbuf_unitized
  const int *_vtIdx_noweld,  // not weld face_sidx, for bvh generation - gltf_model.indexBuffer
  const int *_vtIdx_weld,    // weld face_sidx, for upload cuda memory - face_sidx
  const int *_weld_vidx_all, // weld vidx_distinct_all
  const Vector3 *_fnbuf,     // face normal 
  const Vector4 *_s_edge,    // edge info [v0, v1, t0, t1]
  const Vector4 *_s_emap,    // emap, for lookup generated visibility
  const int *_face_eidx, 
  int _nf,int _dnv, int _nv, int _ne, 
  int* des_ei, float* des_es,
  int* des_vinfo0, 
  int* vis_size)
{
  std::vector<glm::vec3> ls_v;
  std::vector<glm::vec3> ls_vn;
  std::vector<MyVertex> vbuf_unitized;    // vmap
  std::vector<uint32_t> indexBuffer;      // _vtIdx_noweld
  std::vector<uint32_t> face_sidx;        // _vtIdx_weld
  std::vector<int> vidx_distinct_all;     //_weld_vidx_all
  std::vector<glm::vec3> fnbuf;
  std::vector<GEdge> s_edge;
  std::vector<Vector4> tts_edge;
  std::vector<GEdge> s_emap;
  std::vector<Vector4> tts_emap;
  std::vector<int> face_eidx;
  int dnv, nv, ne, nf;  
 
  //copy data to function
  {
    //copy vbuf_unitized
    dnv = _dnv; //should be 1849
    ls_v.resize(dnv);
    ls_vn.resize(dnv); 
    vbuf_unitized.resize(dnv);
    memcpy(ls_v.data(),_ls_v,dnv * sizeof(glm::vec3));
    memcpy(ls_vn.data(),_ls_vn,dnv * sizeof(glm::vec3));

    for(int i=0; i<dnv; i++)
    {
      vbuf_unitized[i].pos = ls_v[i];
      vbuf_unitized[i].normal = ls_vn[i];
    }

    //copy face_sidx
    nf = _nf;
    indexBuffer.resize(nf*3);
    memcpy(indexBuffer.data(), _vtIdx_noweld, nf*3 * sizeof(uint32_t));

    face_sidx.resize(nf*3);
    memcpy(face_sidx.data(),_vtIdx_weld, nf*3 * sizeof(uint32_t));

    //obtain vidx_distinct_all
    nv = _nv;
    vidx_distinct_all.resize(nv);
    memcpy(vidx_distinct_all.data(),_weld_vidx_all, nv* sizeof(int));

    //copy edge info
    ne = _ne;
    tts_edge.resize(ne);
    memcpy(tts_edge.data(),_s_edge,ne * sizeof(Vector4));
    //s_edge.resize(ne);
    for(int i=0;i<tts_edge.size(); i++)
    {
      GEdge aa;
        aa.v0 = int(tts_edge[i].x);
        aa.v1 = int(tts_edge[i].y);
        aa.t0 = int(tts_edge[i].z);
        aa.t1 = int(tts_edge[i].w);
      s_edge.push_back(aa);
    }

    tts_emap.resize(ne);
    memcpy(tts_emap.data(), _s_emap, _ne*sizeof(Vector4));
    for(int i=0; i<tts_emap.size(); i++)
    {
      GEdge aa;
        aa.v0 = int(tts_emap[i].x);
        aa.v1 = int(tts_emap[i].y);
        aa.t0 = int(tts_emap[i].z);
        aa.t1 = int(tts_emap[i].w);
      s_emap.push_back(aa);
    }

    face_eidx.resize(nf*3);
    memcpy(face_eidx.data(),_face_eidx,nf*3 * sizeof(int));

    //copy fnbuf
    fnbuf.resize(nf);
    memcpy(fnbuf.data(), _fnbuf, nf*sizeof(glm::vec3));
  }

  //generate BVH
   bvh0_update(vbuf_unitized.data(), indexBuffer ); 

  VisGen visgen;
	visgen.init_cuda_memory(
		(int)vidx_distinct_all.size(),
		(int*)face_sidx.data(),
		(int*)face_eidx.data(), (int)face_sidx.size() / 3,
		(GEdge*)s_edge.data(), ne
	);
  visgen.update_cuda_memoryX(
    vbuf_unitized.data(), vidx_distinct_all,
    (FLOAT3*)fnbuf.data(), (int)fnbuf.size(), ne
  );

  std::vector<PointI> assembled_dat0; //the generated visibility 
  std::vector<int> assembled_vinfo0;
  int assembled_dat0_buf_size = 0;
  int assembled_vinfo0_buf_size = 0;
  int total_element, n_vertex;

  //evaluate visibility in CUDA
	cal_visibility_rendering_time_ac(
		visgen, vbuf_unitized, vidx_distinct_all, fnbuf,
		ne,
		assembled_dat0, assembled_dat0_buf_size,
		assembled_vinfo0, assembled_vinfo0_buf_size, total_element, n_vertex);



  std::vector<int> ttedge;
  for(int i=0; i<s_edge.size(); i++)
  {
    ttedge.push_back(s_emap[i].v0);
    ttedge.push_back(s_emap[i].t0);
    ttedge.push_back(s_emap[i].v1);
    ttedge.push_back(s_emap[i].t1);
  }

  //copy data to Unity
  int nvinfo = n_vertex+1;    
  std::vector<int> vsize;
    vsize.push_back(total_element);
    vsize.push_back(nvinfo);
  
  std::vector<int> ttei;
  std::vector<float> ttes;
  for(int i=0; i<total_element; i++)
  {
    ttei.push_back(assembled_dat0[i].ei);
    ttes.push_back(assembled_dat0[i].es);
  }

  memcpy(vis_size,vsize.data(),vsize.size() * sizeof(int)); // copy the size of the above three arrays to output
  memcpy(des_ei, ttei.data(), total_element * sizeof(int)); // copy the temproal vec. vis. to output
  memcpy(des_es, ttes.data(), total_element * sizeof(float)); // copy the temproal vec. vis. to output
  memcpy(des_vinfo0,assembled_vinfo0.data(), nvinfo * sizeof(int)); // copy the temproal vec. vis. to output
}


extern "C" __declspec(dllexport) void __stdcall
get_singleV_vecvis_ac(
	const Vector3 * _ls_v, 
  const Vector3 * _ls_vn, // for vbuf_unitized
	const int* _indexbuf, 
  const Vector4* _emap,
	int* _des_ei, float* _des_es,
	int* _des_vinfo0,
	int _nv, int _nf, int _nei, int _dnv, int _nemap,
	int _vidx,
	Vector3 * des_vis, int* n_desvis)
{
	std::vector<glm::vec3> ls_v;
	std::vector<glm::vec3> ls_vn;
	std::vector<int> indexbuf;
	std::vector<Vector4> emap;
	std::vector<int> des_ei;
	std::vector<float> des_es;
	std::vector<int> des_vinfo0;

	int dnv, nv, nf, vidx, nei, nvinfo, nemap;
	  dnv = _dnv;
	  nv = _nv;
	  nf = _nf;
	  nei = _nei;
	  nvinfo = nv + 1;
	  vidx = _vidx;
	  nemap = _nemap;

	//copy data to dll
	ls_v.resize(_dnv);
	  memcpy(ls_v.data(), _ls_v, _dnv * sizeof(glm::vec3));
	ls_vn.resize(_dnv);
	  memcpy(ls_vn.data(), _ls_vn, _dnv * sizeof(glm::vec3));
	indexbuf.resize(nf * 3);
	  memcpy(indexbuf.data(), _indexbuf, nf * 3 * sizeof(int));
	emap.resize(nemap);
	  memcpy(emap.data(), _emap, nemap * sizeof(Vector4));
	des_ei.resize(nei);
	  memcpy(des_ei.data(), _des_ei, nei * sizeof(int));
	des_es.resize(nei);
	  memcpy(des_es.data(), _des_es, nei * sizeof(float));
	des_vinfo0.resize(nvinfo);
	  memcpy(des_vinfo0.data(), _des_vinfo0, nvinfo * sizeof(int));

	int vi = vidx;
  std::vector<glm::vec3> output_vis;
  std::vector<int> output_vissize;
	std::vector<EdgeI> edat;
	for (int i = 0; i < des_ei.size(); i += 2)
	{
		EdgeI aa;
		  aa.ei0 = des_ei[i];//t0.ei;
		  aa.es0 = des_es[i];//t0.es;
		  aa.ei1 = des_ei[i + 1];//t1.ei;
		  aa.es1 = des_es[i + 1];//t1.es;
		edat.push_back(aa);
	}
  uint i0 = des_vinfo0[vidx];
  uint i1 = des_vinfo0[vidx+1];

  glm::vec3 cv = ls_v[ indexbuf[vidx]];
  glm::vec3 ay1 = ls_vn[ indexbuf[vidx]];
  glm::vec3 bx, by, bz;
  get_coordinate_system(ay1, bx, by, bz);
  int i, s;
  for( i=i0; i<i1; i+=2 )
  {
    vec3 u0, u1;
    gete_dll_ac(u0, u1, i/2, cv, bx, bz, ls_v, emap, edat);
    output_vis.push_back(u0);
    output_vis.push_back(u1);
  }
  output_vissize.push_back( output_vis.size());

  memcpy(n_desvis, output_vissize.data(),output_vissize.size() * sizeof(int)); // copy the size of the above three arrays to output
  memcpy(des_vis, output_vis.data(), output_vis.size() * sizeof(glm::vec3)); // copy the temproal vec. vis. to output
}

extern "C" __declspec(dllexport) void __stdcall
gen_tempVecVis_ondemand_ac(
  const Vector3 *_ls_v, 
  const Vector3 *_ls_vn, // for vbuf_unitized
  const int *_vtIdx_noweld,  // not weld face_sidx, for bvh generation - gltf_model.indexBuffer
  const int *_vtIdx_weld,    // weld face_sidx, for upload cuda memory - face_sidx
  const int *_weld_vidx_all, // weld vidx_distinct_all
  const Vector3 *_fnbuf,     // face normal 
  const Vector4 *_s_edge,    // edge info [v0, v1, t0, t1]
  const Vector4 *_s_emap,    // emap, for lookup generated visibility
  const int *_face_eidx, 
  int _nf,int _dnv, int _nv, int _ne, 
  int _first_vi, 
  int* des_ei, float* des_es,
  int* des_vinfo0, 
  int* vis_size)
{
  std::vector<glm::vec3> ls_v;
  std::vector<glm::vec3> ls_vn;
  std::vector<MyVertex> vbuf_unitized;    // vmap
  std::vector<uint32_t> indexBuffer;      // _vtIdx_noweld
  std::vector<uint32_t> face_sidx;        // _vtIdx_weld
  std::vector<int> vidx_distinct_all;     //_weld_vidx_all
  std::vector<glm::vec3> fnbuf;
  std::vector<GEdge> s_edge;
  std::vector<Vector4> tts_edge;
  std::vector<GEdge> s_emap;
  std::vector<Vector4> tts_emap;
  std::vector<int> face_eidx;
  int dnv, nv, ne, nf;  
  int first_vi;
 
  {
    //copy vbuf_unitized
    dnv = _dnv; //should be 1849
    ls_v.resize(dnv);
    ls_vn.resize(dnv); 
    vbuf_unitized.resize(dnv);
    memcpy(ls_v.data(),_ls_v,dnv * sizeof(glm::vec3));
    memcpy(ls_vn.data(),_ls_vn,dnv * sizeof(glm::vec3));

    for(int i=0; i<dnv; i++)
    {
      vbuf_unitized[i].pos = ls_v[i];
      vbuf_unitized[i].normal = ls_vn[i];
    }

    //copy face_sidx
    nf = _nf;
    indexBuffer.resize(nf*3);
    memcpy(indexBuffer.data(), _vtIdx_noweld, nf*3 * sizeof(uint32_t));

    face_sidx.resize(nf*3);
    memcpy(face_sidx.data(),_vtIdx_weld, nf*3 * sizeof(uint32_t));
  
    //obtain vidx_distinct_all
    nv = _nv;
    vidx_distinct_all.resize(nv);
    memcpy(vidx_distinct_all.data(),_weld_vidx_all, nv* sizeof(int));

    //copy edge info
    ne = _ne;
    tts_edge.resize(ne);
    memcpy(tts_edge.data(),_s_edge,ne * sizeof(Vector4));
    for(int i=0;i<tts_edge.size(); i++)
    {
      GEdge aa;
        aa.v0 = int(tts_edge[i].x);
        aa.v1 = int(tts_edge[i].y);
        aa.t0 = int(tts_edge[i].z);
        aa.t1 = int(tts_edge[i].w);
      s_edge.push_back(aa);
    }

    tts_emap.resize(ne);
    memcpy(tts_emap.data(), _s_emap, _ne*sizeof(Vector4));
    for(int i=0; i<tts_emap.size(); i++)
    {
      GEdge aa;
        aa.v0 = int(tts_emap[i].x);
        aa.v1 = int(tts_emap[i].y);
        aa.t0 = int(tts_emap[i].z);
        aa.t1 = int(tts_emap[i].w);
      s_emap.push_back(aa);
    }

    face_eidx.resize(nf*3);
    memcpy(face_eidx.data(),_face_eidx,nf*3 * sizeof(int));

    fnbuf.resize(nf);
    memcpy(fnbuf.data(), _fnbuf, nf*sizeof(glm::vec3));

    first_vi = _first_vi;
  }

  //generate BVH
  bvh0_update(vbuf_unitized.data(), indexBuffer ); 

  VisGen visgen;
	visgen.init_cuda_memory(
		(int)vidx_distinct_all.size(),
		(int*)face_sidx.data(),
		(int*)face_eidx.data(), (int)face_sidx.size() / 3,
		(GEdge*)s_edge.data(), ne
	);
  visgen.update_cuda_memoryX(
    vbuf_unitized.data(), vidx_distinct_all,
    (FLOAT3*)fnbuf.data(), (int)fnbuf.size(), ne
  );


  PointI* assembled_dat0 = NULL;
  int* assembled_vinfo0 = NULL;
  int assembled_dat0_buf_size = 0;
  int assembled_vinfo0_buf_size = 0;
  int total_element, n_vertex;

  if(first_vi!=0)
  {
    // copy original visibility info
    std::vector<int> ass_size;
      ass_size.resize(2);
      memcpy(ass_size.data(), vis_size, 2*sizeof(int));

    std::vector<int> pre_des_ei;
    std::vector<float> pre_des_es;
      pre_des_ei.resize(ass_size[0]);
      pre_des_es.resize(ass_size[0]);
      memcpy(pre_des_ei.data(), des_ei,ass_size[0]*sizeof(int));
      memcpy(pre_des_es.data(), des_es,ass_size[0]*sizeof(float));

    std::vector<PointI> ttassembled_dat0; //the generated visibility 
    std::vector<int> ttassembled_vinfo0;

    for(int j=0; j<pre_des_ei.size(); j++)
    {
      PointI aa;
        aa.ei = pre_des_ei[j];
        aa.es = pre_des_es[j];
      ttassembled_dat0.push_back(aa);
    }
    assembled_dat0 = (PointI*)malloc(pre_des_ei.size()*16*sizeof(PointI));
    memcpy(assembled_dat0, ttassembled_dat0.data(), ttassembled_dat0.size()*sizeof(PointI));

    assembled_vinfo0 = (int*)malloc((nv+1)*16*sizeof(int));
    memcpy(assembled_vinfo0, des_vinfo0, (nv+1)*sizeof(int) );

    assembled_dat0_buf_size = pre_des_ei.size()*sizeof(PointI)*16;
    assembled_vinfo0_buf_size = (nv + 1) * sizeof(int)*16;
  }

  //evaluate visibility in CUDA
  visgen.total_line_sample = 0;
	visgen.update_bvh(bvh0_to_cubvh(), bvh0_node_count());
  {
		int assembled_dat0_size, assembled_vinfo0_size;
		int ipass, npass, total_vertex; 

		if (assembled_vinfo0)
		{
			assembled_dat0_size = assembled_vinfo0[first_vi] * sizeof(PointI);
			assembled_vinfo0_size = (first_vi + 1) * sizeof(int);
			total_element = assembled_vinfo0[first_vi];
		}
		else
		{
			first_vi = 0;
			assembled_dat0_size = 0;
			assembled_vinfo0_size = 0;
			total_element = 0;
		}

		total_vertex = nv - first_vi;
		npass = (total_vertex + (visgen.n_vertex_perpass - 1)) / visgen.n_vertex_perpass;
		n_vertex = first_vi;

		for (ipass = 0; ipass < npass; ipass++)
		{
			int ipass_nv = g_min(total_vertex - ipass * visgen.n_vertex_perpass, visgen.n_vertex_perpass);
			visgen.gen_vis_rendering3(ne, n_vertex, ipass_nv);

			PointI* dat; int dat_size;
			int* vinfo; int vinfo_size;

			dat = visgen.my_dat_buf;
			dat_size = visgen.vis_tightpacked_size * sizeof(PointI);
			vinfo = visgen.my_dat_vi;
			vinfo_size = (ipass_nv + 1) * sizeof(int);

			if (assembled_dat0_buf_size < assembled_dat0_size + dat_size)
			{
				assembled_dat0_buf_size += dat_size * 16;
				assembled_dat0 = (PointI*)realloc(assembled_dat0, assembled_dat0_buf_size);
			}
			if (assembled_vinfo0_buf_size < assembled_vinfo0_size + vinfo_size)
			{
				assembled_vinfo0_buf_size += vinfo_size * 16;
				assembled_vinfo0 = (int*)realloc(assembled_vinfo0, assembled_vinfo0_buf_size);
				assembled_vinfo0[0] = 0;
			}

			memcpy(((BYTE*)assembled_dat0) + assembled_dat0_size, dat, dat_size);
			assembled_dat0_size += dat_size;
			total_element += visgen.my_dat_vi[ipass_nv];

			int i;
			for (i = 1; i <= ipass_nv; i++)
				assembled_vinfo0[n_vertex + i] = assembled_vinfo0[n_vertex] + vinfo[i];
			assembled_vinfo0_size += vinfo_size;
			n_vertex += ipass_nv;

		}
  }

  //copy data to Unity
  int nvinfo = n_vertex+1;    
  std::vector<int> vsize;
    vsize.push_back(total_element);
    vsize.push_back(nvinfo);
  

  std::vector<int> ttei;
  std::vector<float> ttes;
  for(int i=0; i<total_element; i++)
  {
    ttei.push_back(assembled_dat0[i].ei);
    ttes.push_back(assembled_dat0[i].es);
  }

  memcpy(vis_size,vsize.data(),vsize.size() * sizeof(int)); // copy the size of the above three arrays to output
  memcpy(des_ei, ttei.data(), total_element * sizeof(int)); // copy the temproal vec. vis. to output
  memcpy(des_es, ttes.data(), total_element * sizeof(float)); // copy the temproal vec. vis. to output
  memcpy(des_vinfo0,assembled_vinfo0, nvinfo * sizeof(int)); // copy the temproal vec. vis. to output


  if(assembled_vinfo0!=NULL)
    free(assembled_vinfo0);
  if(assembled_dat0!=NULL)
    free(assembled_dat0);
}



//_ls_v and _ls_vn are the vertex position and vertex normal 
//_vtIdx is the triangle vertex index array 
//_dnv is the number of vertex 
//_nf is the number of triangle 
//_weld_vidx_all is weld triangle vertex index array
//_weld_vidx_AllinOne stores the initial vector<vector> vidx_distinct_all in a single array
//_weld_vidx_perSize: each value stores the number of vertex index (for _weld_vidx_AllinOne) which belongs to the same vertex position
extern "C" __declspec(dllexport) void __stdcall
weld_func_ab( 
  const Vector3 *_ls_v,
  const int *_vtIdx_noweld,
  int _dnv, int _nf,
  int *_vtIdx_weld,
  int *_weld_vidx_all,
  int *_weld_vidx_AllinOne,
  int *_weld_vidx_perSize,
  int *_weld_vidx_mapping,
  int *_weldv_size
){
  std::vector<glm::vec3> ls_v;
  std::vector<uint32_t> face_vidx;
  std::vector<int> weld_vidx_mapping;
     

  int dnv, nf;
    dnv = _dnv;
    nf  = _nf;

	ls_v.resize(_dnv);
	  memcpy(ls_v.data(), _ls_v, _dnv * sizeof(glm::vec3));
	face_vidx.resize(nf * 3);
	  memcpy(face_vidx.data(), _vtIdx_noweld, nf * 3 * sizeof(uint32_t));

  weld_vidx_mapping.resize(dnv);

  bool vfound, ifound;
  int i, j, idx;
  uint32_t* v;
  int n;
  int n_vertex;

  std::vector< std::vector<uint32_t> > weld_vidx_all;
  n_vertex = dnv;
  weld_vidx_all.clear();
  
  for( j = 0; j < n_vertex; j++ )
  {
    ifound = false;
    n = (int) face_vidx.size();
    v = face_vidx.data();
    for( i = 0; i < n; i++, v++ )
    {
      if (*v == j)
      {
        ifound = true;
        break;
      }
    }

    if (ifound)
    {
      vfound = false;
      float nd = FLT_MAX;
      for( i = 0; i < weld_vidx_all.size(); i++ )
      {
        float d = length(ls_v[j] - ls_v[weld_vidx_all[i][0]]);
        if (nd > d)
          nd = d;
        if (d < 10 * FLT_EPSILON)
        {
          vfound = true;
          break;
        }
      }

      if( vfound )
      {
        idx = i;
      }
      else
      {
        idx = (int) weld_vidx_all.size();
        weld_vidx_all.push_back( std::vector<uint32_t>() );
      }
      weld_vidx_all[idx].push_back(j);

      n = (int) face_vidx.size();
      v = face_vidx.data();
      for( i=0; i<n; i++, v++ )
      {
        if(*v == j)
          *v = idx;
      }
    }
  }

  std::vector<int> weld_vidx_00;
  std::vector<int> weld_vidx_01;
  std::vector<int> weld_vidx_02;
  std::vector<int> weld_vsize;

  for(i=0; i<weld_vidx_all.size(); i++)
  {
    weld_vidx_00.push_back(weld_vidx_all[i][0]);
    weld_vidx_02.push_back(weld_vidx_all[i].size());

    weld_vidx_mapping[weld_vidx_all[i][0]] = weld_vidx_all[i][0]; //weld intial is the initial


    for(j=0;j<weld_vidx_all[i].size(); j++)
    {
      weld_vidx_01.push_back(weld_vidx_all[i][j]);
      weld_vidx_mapping[weld_vidx_all[i][j]] = i; //weld_vidx_all[i][0];
    }
  }

  std::vector<int> ttsize;
    ttsize.push_back(weld_vidx_00.size());
    ttsize.push_back(weld_vidx_01.size());
    ttsize.push_back(weld_vidx_02.size());

  memcpy(_weld_vidx_all, weld_vidx_00.data(), weld_vidx_00.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_weld_vidx_AllinOne, weld_vidx_01.data(), weld_vidx_01.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_weld_vidx_perSize, weld_vidx_02.data(), weld_vidx_02.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_weldv_size, ttsize.data(), ttsize.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_vtIdx_weld, face_vidx.data(), face_vidx.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_weld_vidx_mapping, weld_vidx_mapping.data(), weld_vidx_mapping.size() * sizeof(int) ); // copy the size of the above three arrays to output

}


extern "C" __declspec(dllexport) void __stdcall
cal_edge_func_ab(
  int *_vtIdx,
  int *_weld_vidx_all,
  int _nv, int _nf, int _weld_nv,
  Vector4 *_s_edge, 
  int* _face_eidx,
  Vector4 *_s_emap,
  int* _nedege)
{
  int n_vertex, nf, weld_nv;
    n_vertex = _nv;
    nf = _nf;
    weld_nv = _weld_nv;

  std::vector<uint32_t> face_vidx;
    face_vidx.resize(nf*3);
    memcpy(face_vidx.data(), _vtIdx, nf * 3 * sizeof(uint32_t));

  std::vector<int> vidx_distinct_all;
    vidx_distinct_all.resize(weld_nv);
    memcpy(vidx_distinct_all.data(), _weld_vidx_all, weld_nv * sizeof(int));

  std::vector<GEdge> s_edge;
  std::vector<int> face_eidx;
  int n_edge;

  if(s_edge.size()!=0)
    s_edge.clear();
  if(face_eidx.size()!=0)
    face_eidx.clear();

  int u[3];
  bool found;
  int i, j, k, ei;
  int s, t;
  int n_face;
  
  n_face = nf;//int( face_vidx.size()/3 );
  s_edge.resize(3 * n_face);
  face_eidx.resize(3 * n_face);
  n_edge = 0;

  std::vector< std::vector<int> > v2e(n_vertex);

  for (i = 0; i < 3 * n_face; i++)
  {
    s_edge[i].t1 = -1;
  }

  for (j = 0; j < n_face; j++)
  {
    u[0] = face_vidx[3 * j + 0];
    u[1] = face_vidx[3 * j + 1];
    u[2] = face_vidx[3 * j + 2];

    for (k = 0; k < 3; k++)
    {
      s = u[k];
      t = u[(k + 1) % 3];
      found = false;

      for (ei = 0; ei < v2e[t].size(); ei++)
      {
        i = v2e[t][ei];
        if (s_edge[i].v1 == s)
        {
          v2e[t].erase(v2e[t].begin() + ei);
          found = true;
          break;
        }
      }

      if (!found)
      {
        s_edge[n_edge].v0 = s;
        s_edge[n_edge].v1 = t;
        s_edge[n_edge].t0 = j;
        face_eidx[3 * j + k] = n_edge;
        v2e[s_edge[n_edge].v0].push_back(n_edge);
        n_edge++;
      }
      else
      {
        s_edge[i].t1 = j;
        face_eidx[3 * j + k] = i;
      }
    }
  }

  bool closed = true;
  for (i = 0; i < n_edge; i++)
  {
    if (s_edge[i].t1 == -1)
    {
      closed = false;
      break;
    }
  }
  s_edge.resize(n_edge);


  std::vector<int> tt;
    tt.push_back(n_edge);

  std::vector<Vector4> ttedge;
  for(i=0;i<n_edge; i++)
  {
    Vector4 aa;
      aa.x =s_edge[i].v0; 
      aa.y =s_edge[i].v1; 
      aa.z =s_edge[i].t0; 
      aa.w =s_edge[i].t1; 
    ttedge.push_back(aa);
  }
  std::vector<Vector4> sedge_lookup;
  for(int i=0; i<s_edge.size(); i++)
  {
    Vector4 aa;
      aa.x = vidx_distinct_all[s_edge[i].v0];
      aa.y = vidx_distinct_all[s_edge[i].v1];
      aa.z = s_edge[i].t0; 
      aa.w = s_edge[i].t1; 
    sedge_lookup.push_back(aa);
  }

  memcpy(_s_edge, ttedge.data(), ttedge.size() * sizeof(Vector4) ); // copy the size of the above three arrays to output
  memcpy(_s_emap, sedge_lookup.data(), sedge_lookup.size() * sizeof(Vector4) ); // copy the size of the above three arrays to output
  memcpy(_face_eidx, face_eidx.data(), face_eidx.size() * sizeof(int) ); // copy the size of the above three arrays to output
  memcpy(_nedege, tt.data(), tt.size() * sizeof(int) ); // copy the size of the above three arrays to output
}

