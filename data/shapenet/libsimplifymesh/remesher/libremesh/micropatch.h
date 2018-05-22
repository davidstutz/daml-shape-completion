#ifndef REMESHER_MICRO_PATCH_HEADER
#define REMESHER_MICRO_PATCH_HEADER

#include <vector>

#include "trianglemesh.h"
#include "defines.h"
#include "vec2.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class calculates a very small patch consisting of a center
 * vertex and neighboring vertices that form a loop around the center
 * vertex. Neighboring vertices need to be specified in order and
 * form a single loop aroung the center vertex.
 *
 * The resulting parametrization to 2D preserves the length of edges
 * and the relative angles between edges incident to the center vertex.
 * Each angle "a_i" results in an angle "alpha * a_i" such that
 *
 *   alpha * sum(a_i) = 2 PI  <==>  alpha = 2 PI / sum(a_i)
 *
 * The position of the center vertex is mapped to the origin.
 */
class MicroPatch
{
  public:
    typedef std::vector<Vec2f> AdjacentVertices2D;
    typedef std::vector<Vec3f> AdjacentVertices3D;

  protected:
    Vec3f center;
    AdjacentVertices3D adj3d;
    AdjacentVertices2D adj2d;

  public:
    /* The center vertex and the vertex neighbors need to be set. */
    void set_center_vertex (Vec3f const& center);
    void append_adjacent_vertex (Vec3f const& ajdacent);

    /* Clears all information to reuse the object for a new patch. */
    void clear (void);

    /* This creates the patch from the vertices. */
    void calculate_patch (void);

    /* After the patch is created, these methods provide information. */
    AdjacentVertices2D const& get_adjacent_vertices (void) const;
    AdjacentVertices2D& get_adjacent_vertices (void);

    /* Creates debug triangle mesh from both, the 2D and the 3D patch. */
    TriangleMeshPtr get_debug_mesh (void) const;
};

/* ---------------------------------------------------------------- */

inline void
MicroPatch::set_center_vertex (Vec3f const& center)
{
  this->center = center;
}

inline void
MicroPatch::append_adjacent_vertex (Vec3f const& adjacent)
{
  this->adj3d.push_back(adjacent);
}

inline MicroPatch::AdjacentVertices2D const&
MicroPatch::get_adjacent_vertices (void) const
{
  return this->adj2d;
}

inline MicroPatch::AdjacentVertices2D&
MicroPatch::get_adjacent_vertices (void)
{
  return this->adj2d;
}

inline void
MicroPatch::clear (void)
{
  this->adj3d.clear();
  this->adj3d.clear();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MICRO_PATCH_HEADER */
