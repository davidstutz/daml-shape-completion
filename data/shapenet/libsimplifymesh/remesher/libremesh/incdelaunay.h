#ifndef REMESHER_INCDELAUNAY_HEADER
#define REMESHER_INCDELAUNAY_HEADER

#include <vector>

#include "defines.h"
#include "vec2.h"
#include "vec3.h"
#include "halfedge.h"

REMESHER_NAMESPACE_BEGIN

/* Internal vertex represenation. */
struct IDVertex
{
  Vec2f pos;
  Vec2f bary;

  IDVertex (Vec2f const& pos, Vec2f const& bary) : pos(pos), bary(bary) {}
};

/* ---------------------------------------------------------------- */

/*
 * Creates an incremental delaunay triangulation for a triangle in 3D.
 * The master triangle is passed to the constructor, samples are
 * inserted with barycentric coordinates w.r.t. the master triangle.
 */
class IncDelaunay : public HalfEdge
{
  public:
    typedef std::vector<IDVertex> IDVertexList;

  private:
    IDVertexList verts;
    float ins_eps;

  private:
    Vec3f get_bary_for_face (HEFace const* face, Vec2f const& pos) const;
    bool is_inside_face (Vec3f const& bary) const;
    void delaunay_insert (HEFace* face);
    void restore_delaunay (HEEdge* edge);

  public:
    IncDelaunay (Vec3f const& v1, Vec3f const& v2, Vec3f const& v3);
    void set_insert_epsilon (float epsilon);

    bool insert_sample (Vec2f const& bary);

    IDVertexList const& get_vertices (void) const;

    void write_debug_mesh (std::string const& filename);
};

/* ---------------------------------------------------------------- */

inline void
IncDelaunay::set_insert_epsilon (float epsilon)
{
  this->ins_eps = epsilon;
}

inline bool
IncDelaunay::is_inside_face (Vec3f const& bary) const
{
  return (bary[0] >= 0.0f && bary[1] >= 0.0f && bary[2] >= 0.0f);
}

inline IncDelaunay::IDVertexList const&
IncDelaunay::get_vertices (void) const
{
  return this->verts;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_INCDELAUNAY_HEADER */
