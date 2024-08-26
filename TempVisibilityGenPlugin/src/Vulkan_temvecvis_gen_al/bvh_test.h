#ifndef BVH_TEST_H
#define BVH_TEST_H


#include <vector>

#include "glm/glm.hpp"
using namespace glm;

#include "g_vector.h"
#include "shaders/my_raytrace_datatypes.h"
#include "psc_aabb.cuh"

struct CloseHitResult
{
  FLOAT3 bc;
  int tidx;
};

void bvh0_update( const MyVertex *vbuf, const std::vector<uint32_t> &ibuf );
cubvh* bvh0_to_cubvh();
int bvh0_node_count();
CloseHitResult* bvh0_raytrace(CloseHitResult *res, FLOAT3 ray_origin, FLOAT3 ray_dir);

void bvh0_to_cubvhVec(std::vector<cubvh> &dat);

#include <bvh/bvh.hpp>
#include <bvh/vector.hpp>
#include <bvh/triangle.hpp>
#include <bvh/ray.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/locally_ordered_clustering_builder.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/primitive_intersectors.hpp>

using Scalar   = float;
using Vector3  = bvh::Vector3<Scalar>;
using Triangle = bvh::Triangle<Scalar>;
using Ray      = bvh::Ray<Scalar>;
using Bvh      = bvh::Bvh<Scalar>;
using Morton = uint32_t;


#endif