#include <iostream>

#include "delaunayflips.h"

REMESHER_NAMESPACE_BEGIN

bool
DelaunayFlips::check_flip_edge (std::size_t v1idx,
    std::size_t v2idx, std::size_t v0idx, std::size_t v3idx,
    std::size_t /* face1 */, std::size_t /* face2 */)
{
  MeshVertexList const& verts = this->mesh->get_vertices();

  /* Check angle: If the sum of the two opposite angles is >= 180, flip! */
  float angle;
  {
    Vec3f e1 = (verts[v1idx] - verts[v0idx]).norm();
    Vec3f e2 = (verts[v2idx] - verts[v0idx]).norm();
    float scalar = e1.scalar(e2);
    scalar = MY_MIN(1.0f, scalar);
    scalar = MY_MAX(-1.0f, scalar);
    angle = std::acos(scalar);

    e1 = (verts[v1idx] - verts[v3idx]).norm();
    e2 = (verts[v2idx] - verts[v3idx]).norm();
    scalar = e1.scalar(e2);
    scalar = MY_MIN(1.0f, scalar);
    scalar = MY_MAX(-1.0f, scalar);
    angle += std::acos(scalar);
  }

  REMESHER_NAN_CHECK(angle)

  /* No need to flip. */
  if (angle < MY_PI || FLOAT_EQ(angle, MY_PI))
    return false;

  /* Check if flipping would destroy topology. */
  if (this->vinfo->is_mesh_edge(v0idx, v3idx))
  {
    //std::cout << "Flip denied! " << v0idx << " -> " << v3idx << std::endl;
    return false;
  }

  return true;
}

REMESHER_NAMESPACE_END
