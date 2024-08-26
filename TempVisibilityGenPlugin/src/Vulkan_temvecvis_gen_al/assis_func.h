#ifndef ASSIS_FUNC_H
#define ASSIS_FUNC_H

#include <vector>
#include "glm/glm.hpp"
#include "shaders/my_raytrace_datatypes.h"

void save_MyVertexArr(char *spath, std::vector<MyVertex> &vbuf, int nv);
void save_intArr(char *spath, std::vector<int> &ibuf, int nbuf, int npl);
void save_vec3Arr(char *spath, std::vector<glm::vec3> &arr, int narr);
void save_uintArr(char *spath, std::vector<uint32_t> &ibuf, int nbuf, int npl);


#endif // !ASSIS_FUNC_H
