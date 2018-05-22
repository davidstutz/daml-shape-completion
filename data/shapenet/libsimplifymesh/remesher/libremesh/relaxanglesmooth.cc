#include "relaxanglesmooth.h"

REMESHER_NAMESPACE_BEGIN

Vec2f
AngleSmoothing::micropatch_relocate (MicroDelaunay const& md)
{
  MicroPatch::AdjacentVertices2D const& verts = md.get_adjacent_vertices();
  Vec2f center(0.0f, 0.0f);

  Vec2f newpos(0.0f, 0.0f);
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    std::size_t im1 = (i + verts.size() - 1) % verts.size();
    std::size_t ip1 = (i + 1) % verts.size();

    Vec2f const& vm1 = verts[im1];
    Vec2f const& vi = verts[i];
    Vec2f const& vp1 = verts[ip1];

    Vec2f pos = (vm1.norm() + vp1.norm()) * vi.length() / 2.0f;
    newpos += pos;
  }

  newpos /= (float)verts.size();
  return newpos;
}

REMESHER_NAMESPACE_END
