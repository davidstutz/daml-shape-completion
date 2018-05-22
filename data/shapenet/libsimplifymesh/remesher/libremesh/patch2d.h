#ifndef REMESHER_PATCH_2D_HEADER
#define REMESHER_PATCH_2D_HEADER

#include <vector>

#include "defines.h"
#include "refptr.h"
#include "patch3d.h"
#include "vertexinfo.h"
#include "trianglemesh.h"

#define CREATE_DEBUG_3D_MESHES 0
#define CREATE_DEBUG_2D_MESHES 0
#define DEBUG_MESH_3D_PATCH_FILE "../debug/patch3d.off"
#define DEBUG_MESH_2D_PATCH_FILE "../debug/patch2d.off"

REMESHER_NAMESPACE_BEGIN

/*
 * This class creates a 2D parametrization of a set of connected
 * faces (3D patch). The faces are expected to be isomorphic to a disc.
 *
 * The boundary of the triangulation is embetted as polygon whose
 * vertices lie on the unit circle. The distances between the vertices
 * on the circle are proportional to those on the reference mesh.
 * A linear system Ax = b with with N * N matrix elements, where
 * N refers to the amount of non-boundary vertices, is solved in
 * order to aquire the new 2D coordinates.
 *
 * In weights used in this implementation are those described in
 *
 *   Michael S. Floater
 *   Mean Value Coordinates
 *
 * The technique produces a conformal mapping, thus minimizing
 * angular distortion, and guarantees a valid parametrization.
 */

class Patch2d;
typedef RefPtr<Patch2d> Patch2dPtr;

class Patch2d
{
  public:
    typedef std::vector<std::size_t> PatchBoundary;
    typedef std::vector<std::size_t> IndexMapping;
    typedef IndexMapping VertexMapping;
    typedef IndexMapping FaceMapping;

  private:
    VertexMapping vmap;
    FaceMapping fmap;
    TriangleMeshPtr mesh3d;
    TriangleMeshPtr mesh2d;

  protected:
    Patch2d (void);

    void create_patch (Patch3dPtr patch);
    void create_patch (TriangleMeshPtr mesh);

    void create_index_mappings (Patch3dPtr patch);
    void create_identity_index_mappings (void);
    void create_patch_mesh (void);
    void create_parametrization (void);

    void create_patch_boundary (PatchBoundary& b, VertexInfoListPtr vinfo);
    void create_embedding (PatchBoundary& b, VertexInfoListPtr vinfo);

    /* Intern logarithmic lookup. */
    std::size_t lookup_intern (IndexMapping const& map, std::size_t i,
        std::size_t a, std::size_t b) const;

    /* Reordering of bouondary vertices. Vertices are reordered in a way that
     * v_1 to v_k are inner vertices and v_k+1, v_n are border vertices for
     * vertices: v_1, ..., v_k, v_k+1, ..., v_n. */
    void reorder_boundary_vertices (VertexInfoListPtr vinfo);

    /* Lookup of vertices in logarithmic time. This works ONLY PRIOR
     * paramatrization, because border vertices are reordered. */
    std::size_t lookup_vertex (std::size_t index) const;

  public:
    /* Parameterize a part of the mesh given by a set of faces. */
    static Patch2dPtr create (Patch3dPtr patch);

    /* Parameterize a whole mesh. Requires isomorphy to a disc. */
    static Patch2dPtr create (TriangleMeshPtr mesh);

    /* Returns the parameterized patch as triangle mesh. */
    TriangleMeshPtr get_mesh2d (void) const;

    /* Lookup of faces in logarithmic time. */
    std::size_t lookup_face (std::size_t index) const;

    VertexMapping const& get_vertex_mapping (void) const;
    FaceMapping const& get_face_mapping (void) const;

    std::size_t get_memory_usage (void) const;
};

/* ---------------------------------------------------------------- */

inline
Patch2d::Patch2d (void)
{
}

inline Patch2dPtr
Patch2d::create (Patch3dPtr patch)
{
  Patch2dPtr ret(new Patch2d);
  ret->create_patch(patch);
  return ret;
}

inline Patch2dPtr
Patch2d::create (TriangleMeshPtr mesh)
{
  Patch2dPtr ret(new Patch2d);
  ret->create_patch(mesh);
  return ret;
}

inline std::size_t
Patch2d::lookup_vertex (std::size_t index) const
{
  if (this->vmap.empty())
    return MAX_SIZE_T;
  return this->lookup_intern(this->vmap, index, 0, this->vmap.size() - 1);
}

inline std::size_t
Patch2d::lookup_face (std::size_t index) const
{
  if (this->fmap.empty())
    return MAX_SIZE_T;
  return this->lookup_intern(this->fmap, index, 0, this->fmap.size() - 1);
}

inline Patch2d::VertexMapping const&
Patch2d::get_vertex_mapping (void) const
{
  return this->vmap;
}

inline Patch2d::FaceMapping const&
Patch2d::get_face_mapping (void) const
{
  return this->fmap;
}

inline TriangleMeshPtr
Patch2d::get_mesh2d (void) const
{
  return this->mesh2d;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PATCH_2D_HEADER */
