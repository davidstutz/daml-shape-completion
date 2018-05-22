#ifndef REMESHER_MESH_SLICING_HEADER
#define REMESHER_MESH_SLICING_HEADER

#include <map>

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "plane.h"

REMESHER_NAMESPACE_BEGIN

struct MeshSlicingConf
{
  float vertex_snapping;
  float edge_snapping;
};

/* ---------------------------------------------------------------- */

class MeshSlicing
{
  protected:
    /* Data structures for immediate slicing. */
    typedef std::map<std::pair<std::size_t, std::size_t>, std::size_t> NewV;

    /* Data structures for lazy slicing. */
    typedef std::vector<size_t> EdgeVertices;
    typedef std::map<std::pair<std::size_t, std::size_t>,
        EdgeVertices > EdgeMap;

  protected:
    MeshSlicingConf conf;
    TriangleMeshPtr mesh;

    /* Data structure for immediate slicing. */
    NewV new_verts;

    /* Data structure for lazy slicing. */
    EdgeMap edge_info;

    /* This methods cuts a single edge against the cut plane. The
     * resulting vertex is returned for immediate triangulation.
     * The returned vertex is unique for the edge (and was possibly
     * created at a previous intersecion operation). */
    std::size_t get_intersection (std::size_t v1idx,
        std::size_t v2idx, Plane3f const& plane);

    /* This method lazy-cuts a single edge against the cut plane.
     * The new vertex is added to the mesh only if a similar
     * vertex for that edge is not present. */
    void process_intersection (std::size_t v1idx,
        std::size_t v2idx, Plane3f const& plane);

  public:
    MeshSlicing (void);

    void set_config (MeshSlicingConf const& config);
    void set_mesh (TriangleMeshPtr mesh);

    /* Slices the mesh with the given plane and immediately
     * creates a new triangulation for each triangle. This results
     * in bad triangles if there are a lot of skinny triangles. */
    void slice (Plane3f const& plane);

    /* Lazy slicing only calculates intersections with edges but
     * does not triangulate immediately. Only new vertices are
     * created and added to the mesh. */
    void lazy_slice (Plane3f const& plane);

    /* After lazy slicing, each triangle in the mesh must be
     * triangulated using the new intersection vertices. */
    void triangulate (void);
};

/* ---------------------------------------------------------------- */

inline
MeshSlicing::MeshSlicing (void)
{
}

inline void
MeshSlicing::set_config (MeshSlicingConf const& config)
{
  this->conf = config;
}

inline void
MeshSlicing::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_SLICING_HEADER */
