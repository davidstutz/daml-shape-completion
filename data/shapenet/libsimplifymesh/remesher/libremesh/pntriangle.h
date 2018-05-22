#ifndef REMESHER_PN_TRIANGLE_HEADER
#define REMESHER_PN_TRIANGLE_HEADER

#include "defines.h"
#include "trianglemesh.h"
#include "featureedges.h"
#include "vec2.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

class PNTriangle
{
  private:
    TriangleMeshPtr mesh;
    FeatureEdgesPtr features;
    VertexInfoListPtr vinfo;

    Vec3f v[3];
    Vec3f n[3];

  public:
    PNTriangle (void);

    /* The mesh is required if init(face_id) is used. */
    void set_mesh (TriangleMeshPtr mesh);
    /* Surface no considered smooth on feature edges if set. */
    void set_features (FeatureEdgesPtr features);
    /* Vertex info is required if features are set. */
    void set_vertex_info (VertexInfoListPtr vinfo);

    /* Requires valid mesh with set_mesh(). */
    void init (std::size_t face_id);
    /* No requirements. */
    void init (Vec3f const& v1, Vec3f const& v2, Vec3f const& v3,
        Vec3f const& n1, Vec3f const& n2, Vec3f const& n3);

    /* Returns PN position given by barycentric coordinate.
     * Requires that one of the init functions has been called. */
    Vec3f get_position (Vec3f const& bary);
    Vec3f get_position (Vec2f const& bary);
};

/* ---------------------------------------------------------------- */

inline
PNTriangle::PNTriangle (void)
{
}

inline void
PNTriangle::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
PNTriangle::set_features (FeatureEdgesPtr features)
{
  this->features = features;
}

inline void
PNTriangle::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
PNTriangle::init (Vec3f const& v1, Vec3f const& v2, Vec3f const& v3,
    Vec3f const& n1, Vec3f const& n2, Vec3f const& n3)
{
  this->v[0] = v1; this->v[1] = v2; this->v[2] = v3;
  this->n[0] = n1; this->n[1] = n2; this->n[2] = n3;
}

inline Vec3f
PNTriangle::get_position (Vec2f const& bary)
{
  return this->get_position(Vec3f(bary[0], bary[1], 1.0f - bary[0] - bary[1]));
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PN_TRIANGLE_HEADER */
