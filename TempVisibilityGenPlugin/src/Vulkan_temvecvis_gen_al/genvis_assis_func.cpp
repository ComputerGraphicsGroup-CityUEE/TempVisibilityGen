#include "genvis_assis_func.h"
#include <Windows.h>
#include "g_common.h"
#include "g_pfm.h"
#include "psc_aabb.cuh"


void update_visibility_ac(VisGen& visgen, std::vector<MyVertex>& vbuf_unitized,
	std::vector<int>& vidx_distinct_all, std::vector<glm::vec3>& fnbuf, int n_edge,
	std::vector<PointI>& ttassembled_dat0, int& assembled_dat0_buf_size,
	std::vector<int>& ttassembled_vinfo0, int& assembled_vinfo0_buf_size,
	int first_vi, int& total_element, int& n_vertex)
{
	visgen.update_cuda_memoryX(vbuf_unitized.data(), vidx_distinct_all, (FLOAT3*)fnbuf.data(), (int)fnbuf.size(), n_edge);
	visgen.update_bvh(bvh0_to_cubvh(), bvh0_node_count());

	int assembled_dat0_size, assembled_vinfo0_size;
	int ipass, npass, total_vertex;// total_element, n_vertex;

	PointI* assembled_dat0; //= &ttassembled_dat0[0];
	int* assembled_vinfo0; //= &ttassembled_vinfo0[0];

	if (ttassembled_dat0.size() != 0)
		assembled_dat0 = &ttassembled_dat0[0];
	else
		assembled_dat0 = NULL;
	if (ttassembled_vinfo0.size() != 0)
		assembled_vinfo0 = &ttassembled_vinfo0[0];
	else
		assembled_vinfo0 = NULL;

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

	total_vertex = vidx_distinct_all.size() - first_vi;
	npass = (total_vertex + (visgen.n_vertex_perpass - 1)) / visgen.n_vertex_perpass;
	n_vertex = first_vi;

	for (ipass = 0; ipass < npass; ipass++)
	{
		int ipass_nv = g_min(total_vertex - ipass * visgen.n_vertex_perpass, visgen.n_vertex_perpass);
		visgen.gen_vis_rendering3(n_edge, n_vertex, ipass_nv);

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

	ttassembled_dat0.resize(total_element);
	ttassembled_vinfo0.resize(n_vertex + 1);
	memcpy(ttassembled_dat0.data(), assembled_dat0, total_element * sizeof(PointI));
	memcpy(ttassembled_vinfo0.data(), assembled_vinfo0, (n_vertex + 1) * sizeof(int));
}


void cal_visibility_rendering_time_ac(
	VisGen& visgen, std::vector<MyVertex>& vbuf_unitized,
	std::vector<int>& vidx_distinct_all, std::vector<glm::vec3>& fnbuf,
	int n_edge,
	std::vector<PointI>& assembled_dat0, int& assembled_dat0_buf_size,
	std::vector<int>& assembled_vinfo0, int& assembled_vinfo0_buf_size,
	int& total_element, int& n_vertex
)
{
	visgen.total_line_sample = 0;
	int first_vi = 0;
	update_visibility_ac(visgen, vbuf_unitized, vidx_distinct_all, fnbuf, n_edge,
		assembled_dat0, assembled_dat0_buf_size, assembled_vinfo0, assembled_vinfo0_buf_size,
		first_vi, total_element, n_vertex);
}


void get_coordinate_system(glm::vec3 ay, glm::vec3& bx, glm::vec3& by, glm::vec3& bz)
{
	by = normalize(ay);
	bx = glm::vec3(1, 0, 0);
	glm::vec3 ref = mix(glm::vec3(by.z, 0.f, -by.x), bx, abs(by.y));
	bz = normalize(cross(ref, by));
	bx = normalize(cross(by, bz));
}

void gete_dll_ac(glm::vec3& v0, glm::vec3& v1, uint ei, glm::vec3 av, glm::vec3 bx, glm::vec3 bz,
	std::vector<glm::vec3>& vpos, std::vector<Vector4>& emap, std::vector<EdgeI>& edat)
{
	int ei0 = edat[ei].ei0;
	if (ei0 != -1)
	{
		glm::vec3 ev0 = vpos[int(emap[ei0].x)].xyz;
		glm::vec3 ev1 = vpos[int(emap[ei0].y)].xyz;
		v0 = edat[ei].es0 * (ev1 - ev0) + (ev1 + ev0) / 2.f;
	}
	else
	{
		v0 = av + (cos(edat[ei].es0) * bz + sin(edat[ei].es0) * bx) * 2.f;
	}

	int ei1 = edat[ei].ei1;
	if (ei1 != -1)
	{
		glm::vec3 ev0 = vpos[int(emap[ei1].x)].xyz;
		glm::vec3 ev1 = vpos[int(emap[ei1].y)].xyz;
		v1 = edat[ei].es1 * (ev1 - ev0) + (ev1 + ev0) / 2.f;
	}
	else
	{
		v1 = av + (cos(edat[ei].es1) * bz + sin(edat[ei].es1) * bx) * 2.f;
	}
}

///////////////////////////////////////////////////////////////////


void save_intArr(char *spath, std::vector<int> &ibuf, int nbuf, int npl)
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

  GPath gp = parse_spath(spath);
  char s0[256];
  sprintf(s0, "%s.pf1", gp.fname);

  int w = 512;
  int h = nbuf/w +1;
  GPf1 des;
  des.load(w,h);
	for(i=0; i<w*h; i++)
  {
    if(i<nbuf)
      des.fm[i] = ibuf[i];
  }
  des.save(s0);
}

void save_floatArr(char *spath, std::vector<float> &ibuf, int nbuf, int npl)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nbuf; i++)
	{
		fprintf(f0, "%f\t",ibuf[i]);
		if((i+1)%npl==0)
			fprintf(f0, "\n");
	}

	fclose(f0);

  GPath gp = parse_spath(spath);
  char s0[256];
  sprintf(s0, "%s.pf1", gp.fname);

  int w = 512;
  int h = nbuf/w +1;
  GPf1 des;
  des.load(w,h);
	for(i=0; i<w*h; i++)
  {
    if(i<nbuf)
      des.fm[i] = ibuf[i];
  }
  des.save(s0);

}

void save_uintArr(char *spath, std::vector<uint32_t> &ibuf, int nbuf, int npl)
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

void save_vec3Arr(char *spath, std::vector<glm::vec3> &arr, int narr)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<narr; i++)
		fprintf(f0, "%f\t%f\t%f\t\n",arr[i].x,arr[i].y,arr[i].z);

	fclose(f0);
}

void save_vec3ArrPFM(char *spath, std::vector<glm::vec3> &arr, int narr)
{
  int w = 512;
  int h = narr/w +1;
  printf("w, h, w*h, narr %d %d %d %d\n", w, h, w*h, narr);
  GPfm des;
  des.load(w,h);
	int i;
	for(i=0; i<w*h; i++)
  {
    if(i<narr)
      des.fm[i] = FLOAT3(arr[i].x,arr[i].y,arr[i].z);
  }
  des.save(spath);
}

void save_Vector4Arr(char *spath, std::vector<Vector4> &arr, int narr)
{
	int i;  
  FILE *f0 = fopen(spath, "wt");
	for(i=0; i<narr; i++)
		fprintf(f0, "%d\t%d\t%d\t%d\t\n",int(arr[i].x),int(arr[i].y),int(arr[i].z),int(arr[i].w) );
	fclose(f0);

  GPath gp = parse_spath(spath);
  char s0[256];
  sprintf(s0, "%s_single.txt", gp.fname);
	FILE *f1 = fopen(s0, "wt");
	for(i=0; i<narr; i++)
  {
		fprintf(f1, "%d\n",int(arr[i].x));
		fprintf(f1, "%d\n",int(arr[i].y));
		fprintf(f1, "%d\n",int(arr[i].z));
		fprintf(f1, "%d\n",int(arr[i].w));
  }
	fclose(f1);

}

void save_MyVertexArr(const char *spath, std::vector<MyVertex> &vbuf, int nv)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nv; i++)
	{
		fprintf(f0, "%f\t%f\t%f\t\n",vbuf[i].pos.x,vbuf[i].pos.y,vbuf[i].pos.z);
		fprintf(f0, "%f\t%f\t%f\t\n",vbuf[i].normal.x,vbuf[i].normal.y,vbuf[i].normal.z);
	}

	fclose(f0);
}

void save_PointIArr(char *spath, std::vector<PointI> &pibuf, int npi)
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

void save_vidx_distinct(char *spath,  std::vector< std::vector<uint32_t> > &vidx_distinct_all)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<vidx_distinct_all.size(); i++)
  {
    for(int j=0; j<vidx_distinct_all[i].size(); j++)
      fprintf(f0, "%d\t", vidx_distinct_all[i][j]);
    fprintf(f0, "\n");
  }

	fclose(f0);
}

void save_GedgeArr(char *spath, std::vector<GEdge> &arr, int narr)
{
	int i;  
  FILE *f0 = fopen(spath, "wt");
	for(i=0; i<narr; i++)
		fprintf(f0, "%d\t%d\t%d\t%d\t\n",arr[i].v0,arr[i].t0,arr[i].v1,arr[i].t1 );
	fclose(f0);

  GPath gp = parse_spath(spath);
  char s0[256];
  sprintf(s0, "%s_single.txt", gp.fname);
  FILE *f1 = fopen(spath, "wt");
	for(i=0; i<narr; i++)
  {
		fprintf(f1, "%d\n",arr[i].v0);
		fprintf(f1, "%d\n",arr[i].v1);
		fprintf(f1, "%d\n",arr[i].t0);
		fprintf(f1, "%d\n",arr[i].t1);
  }
	fclose(f1);


}
