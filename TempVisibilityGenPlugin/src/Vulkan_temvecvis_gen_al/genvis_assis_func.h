#ifndef GENVIS_ASSIS_FUNC_H
#define GENVIS_ASSIS_FUNC_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"

#include "shaders/my_raytrace_datatypes.h"
#include "cal_visibility_cuda.h"
#include "Qtype.h"

void update_visibility_ac(VisGen& visgen, std::vector<MyVertex>& vbuf_unitized,
	std::vector<int>& vidx_distinct_all, std::vector<glm::vec3>& fnbuf, int n_edge,
	std::vector<PointI>& ttassembled_dat0, int& assembled_dat0_buf_size,
	std::vector<int>& ttassembled_vinfo0, int& assembled_vinfo0_buf_size,
	int first_vi, int& total_element, int& n_vertex);

void cal_visibility_rendering_time_ac(
	VisGen& visgen, std::vector<MyVertex>& vbuf_unitized,
	std::vector<int>& vidx_distinct_all, std::vector<glm::vec3>& fnbuf, int n_edge,
	std::vector<PointI>& assembled_dat0, int& assembled_dat0_buf_size,
	std::vector<int>& assembled_vinfo0, int& assembled_vinfo0_buf_size,
	int& total_element, int& n_vertex );

void get_coordinate_system(glm::vec3 ay, glm::vec3& bx, glm::vec3& by, glm::vec3& bz);

void gete_dll_ac(glm::vec3& v0, glm::vec3& v1, uint ei, glm::vec3 av, glm::vec3 bx, glm::vec3 bz,
	std::vector<glm::vec3>& vpos, std::vector<Vector4>& emap, std::vector<EdgeI>& edat);

///////////////////////////////////////////////////////////////////
void save_intArr(char *spath, std::vector<int> &ibuf, int nbuf, int npl);
void save_floatArr(char *spath, std::vector<float> &ibuf, int nbuf, int npl);
void save_uintArr(char *spath, std::vector<uint32_t> &ibuf, int nbuf, int npl);
void save_vec3Arr(char *spath, std::vector<glm::vec3> &arr, int narr);
void save_vec3ArrPFM(char *spath, std::vector<glm::vec3> &arr, int narr);
void save_Vector4Arr(char *spath, std::vector<Vector4> &arr, int narr);
void save_MyVertexArr(const char *spath, std::vector<MyVertex> &vbuf, int nv);
void save_PointIArr(char *spath, std::vector<PointI> &pibuf, int npi);
void save_vidx_distinct(char *spath,  std::vector< std::vector<uint32_t> > &vidx_distinct_all);
void save_GedgeArr(char *spath, std::vector<GEdge> &arr, int narr);

#endif // !GENVIS_ASSIS_FUNC_H

