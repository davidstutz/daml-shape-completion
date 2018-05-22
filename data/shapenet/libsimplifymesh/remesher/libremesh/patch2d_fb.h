#ifndef REMESHER_PATCH_2D_FB_HEADER
#define REMESHER_PATCH_2D_FB_HEADER

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
 * The parametrization technique that is used here is called linABF,
 * a free boundary parametrization that requires to solve one linear
 * system of equations.
 */

class Patch2dFB;
typedef RefPtr<Patch2dFB> Patch2dFBPtr;

// Temp
typedef Patch2dFB Patch2d;
typedef Patch2dFBPtr Patch2dPtr;

class Patch2dFB
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
	
    std::vector<double> optimalAngles;
    std::vector<double> alphas;


  protected:
    Patch2dFB (void);

    void create_patch (Patch3dPtr patch);
    void create_patch (TriangleMeshPtr mesh);

    void create_index_mappings (Patch3dPtr patch);
    void create_identity_index_mappings (void);
    void create_patch_mesh (void);

    void create_patch_boundary (PatchBoundary& b, VertexInfoListPtr vinfo);
    void create_parametrization (void);
    void create_embedding (VertexInfoListPtr vinfo);
    void linearAnglesToMesh();
    
    void compute_optimalAngles(VertexInfoListPtr vinfo);
    float getAngleSum3d(VertexInfoListPtr vinfo, std::size_t vId);
    float getAngle3D(std::size_t t, std::size_t k);
    void anglesToMesh(VertexInfoListPtr vinfo, std::size_t v1, std::size_t v2);
	

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
    static Patch2dFBPtr create (Patch3dPtr patch);

    /* Parameterize a whole mesh. Requires isomorphy to a disc. */
    static Patch2dFBPtr create (TriangleMeshPtr mesh);

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
Patch2dFB::Patch2dFB (void)
{
}

inline Patch2dFBPtr
Patch2dFB::create (Patch3dPtr patch)
{
  Patch2dFBPtr ret(new Patch2dFB);
  ret->create_patch(patch);
  return ret;
}

inline Patch2dFBPtr
Patch2dFB::create (TriangleMeshPtr mesh)
{
  Patch2dFBPtr ret(new Patch2dFB);
  ret->create_patch(mesh);
  return ret;
}

inline std::size_t
Patch2dFB::lookup_vertex (std::size_t index) const
{
  if (this->vmap.empty())
    return MAX_SIZE_T;
  return this->lookup_intern(this->vmap, index, 0, this->vmap.size() - 1);
}

inline std::size_t
Patch2dFB::lookup_face (std::size_t index) const
{
  if (this->fmap.empty())
    return MAX_SIZE_T;
  return this->lookup_intern(this->fmap, index, 0, this->fmap.size() - 1);
}

inline Patch2dFB::VertexMapping const&
Patch2dFB::get_vertex_mapping (void) const
{
  return this->vmap;
}

inline Patch2dFB::FaceMapping const&
Patch2dFB::get_face_mapping (void) const
{
  return this->fmap;
}

inline TriangleMeshPtr
Patch2dFB::get_mesh2d (void) const
{
  return this->mesh2d;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PATCH_2D_FB_HEADER */
