#include <iostream>

#include "exception.h"
#include "patch3d.h"
#include "relocation.h"

#define MAX_SEARCH_ITERATIONS 100

REMESHER_NAMESPACE_BEGIN

VertexRef
Relocation::relocate (VertexRef const& v1, VertexRef const& v2,
    VertexRef const& v3, Vec2f const& bary_coords)
{
  /* Get a patch for the reference faces from the cache. */
  Patch2dPtr patch2d = this->get_patch(v1.face, v2.face, v3.face);
  TriangleMeshPtr mesh2d = patch2d->get_mesh2d();

  /* Now find the new triangle. Steps are as follows:
   * - Calculate the location of the point, that is
   *     p0 = b1 p1 + b2 p2 + b3 p3
   *   where pn is the 2D-coordinate of the vertex reference
   * - Pick a starting triangle, f = v1.face
   * - Iterate:
   *   - Calculate barycentric coordinates of p0 regarding f
   *   - If all coordiantes are positive, we found the triangle
   *   - Otherwise advance f in the direction of the point
   */

  std::size_t v1_idx = patch2d->lookup_face(v1.face);
  std::size_t v2_idx = patch2d->lookup_face(v2.face);
  std::size_t v3_idx = patch2d->lookup_face(v3.face);
  Triangle2f tri1(mesh2d, v1_idx);
  Triangle2f tri2(mesh2d, v2_idx);
  Triangle2f tri3(mesh2d, v3_idx);

  /* Calculate the 2D point p1, p2, p3. */
  Vec2f p1 = tri1.get_point_for_bary(v1.bary);
  Vec2f p2 = tri2.get_point_for_bary(v2.bary);
  Vec2f p3 = tri3.get_point_for_bary(v3.bary);
  Triangle2f big_tri(p1, p2, p3);
  Vec2f p0 = big_tri.get_point_for_bary(bary_coords);

  #if 0
  std::cout << "Lookup results: " << v1_idx << " " << v2_idx << " " << v3_idx << std::endl;
  std::cout << "Triangle A: " << tri1.a << " " << tri1.b << " " << tri1.c << std::endl;
  std::cout << "Bary: " << v1.bary << std::endl;

  std::cout << "Tri1 Point 1: " << (tri1.a * v1.bary[0] + tri1.b * v1.bary[1] + tri1.c * (1.0f - v1.bary[0] - v1.bary[1])) << std::endl;
  std::cout << "Tri1 Point 2: " << p1 << std::endl;

  std::cout << p1[0] << " " << p1[1] << std::endl;
  std::cout << p2[0] << " " << p2[1] << std::endl;
  std::cout << p3[0] << " " << p3[1] << std::endl;
  std::cout << p0[0] << " " << p0[1] << std::endl;
  #endif

  /* Pick starting face and search the new face. */
  std::size_t f = v1_idx;
  std::size_t iterations = 0;
  bool found = false;
  Vec3f rbary;

  while (!found && iterations < MAX_SEARCH_ITERATIONS)
  {
    /* Calculate barycentric coordinates for current face
     * and check if point is inside. */
    Triangle2f tri(mesh2d, f);
    rbary = tri.get_bary_coords(p0);
    //std::cout << "Barycentric coords for patch face " << f
    //    << ": " << rbary << std::endl;
    if (rbary[0] >= 0.0f && rbary[1] >= 0.0f && rbary[2] >= 0.0f)
    {
      found = true;
      break;
    }
    iterations += 1;

    /* Barycentric coordinates reveal that the point is outside
     * the current triangle. Advance to the next triangle. */
    float smallest_bary;
    std::size_t smallest_idx;
    for (int i = 0; i < 3; ++i)
    {
      if (i == 0 || rbary[i] < smallest_bary)
      {
        smallest_bary = rbary[i];
        smallest_idx = i;
      }
    }

    /* Since no vertex information is available for the patch, the
     * faces of the reference mesh are used to find the next triangle.
     * This corresponding triangle in 2D needs to be looked up. */
    std::size_t ref_face = patch2d->get_face_mapping()[f];
    std::size_t new_face = this->find_opposite_face(ref_face, smallest_idx);
    f = patch2d->lookup_face(new_face);
    if (f == MAX_SIZE_T)
      throw Exception("Face search reached patch boundary!");

    #if 0
    std::cout << "Reference face: " << ref_face
      << " New face: " << new_face
      << " for vertex " << this->rmesh->get_faces()[ref_face * 3 + smallest_idx]
      << " new patch face: " << f << std::endl;
    #endif
  }

  if (iterations >= MAX_SEARCH_ITERATIONS)
    throw Exception("Triangle search failed!");

  VertexRef ret;
  ret.face = patch2d->get_face_mapping()[f];
  ret.bary = Vec2f(rbary[0], rbary[1]);
  return ret;
}

/* ---------------------------------------------------------------- */

Patch2dPtr
Relocation::get_patch (std::size_t f1, std::size_t f2, std::size_t f3)
{
#if REMESHER_PATCH_CACHING
  /* If patch caching is enabled, try to lookup the patch in the cache. */
  {
    Patch2dPtr patch2d = this->cache->lookup_patch(f1, f2, f3);
    if (patch2d.get() != 0)
      return patch2d;
  }
#endif

  /*
   * No patch has been found in the cache, create a new one.
   * First, create the 3D patch fron the three reference triangles.
   */
  Patch3dTrianglePtr patch3d;
  try
  {
    patch3d = Patch3dTriangle::create(this->rmesh, this->vinfo);
    patch3d->create_patch(f1, f2, f3);
  }
  catch (Exception& e)
  {
    throw Exception("Error creating 3D patch: " + e);
  }

  /* Flatten the 3D patch to 2D. */
  Patch2dPtr patch2d;
  try
  {
    patch2d = Patch2d::create((Patch3dPtr)patch3d);
  }
  catch (Exception& e)
  {
    throw Exception("Error creating 2D patch: " + e);
  }

#if REMESHER_PATCH_CACHING
  /* If patch caching is enabled, cache the patch. */
  this->cache->cache_patch(patch2d);
#endif

  return patch2d;
}

/* ---------------------------------------------------------------- */

std::size_t
Relocation::find_opposite_face (std::size_t face, std::size_t fvidx)
{
  MeshFaceList const& faces = this->rmesh->get_faces();

  std::size_t prev_fvidx = (fvidx + 2) % 3;
  std::size_t vidx = faces[face * 3 + prev_fvidx];
  VertexInfo::FaceList const& adj_faces = this->vinfo[vidx].adj_faces;

  for (std::size_t i = 0; i < adj_faces.size(); ++i)
    if (adj_faces[i] == face)
      return adj_faces[(i + 1) % adj_faces.size()];

  throw Exception("Cannot find opposite face");
}

REMESHER_NAMESPACE_END
