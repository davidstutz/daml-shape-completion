#ifndef REMESHER_MESH_DECIMATION_HEADER
#define REMESHER_MESH_DECIMATION_HEADER

#include <vector>

#include "defines.h"
#include "trianglemesh.h"
#include "featureedges.h"
#include "meshcleanup.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

/* The mesh decimation algorithm takes as input a mesh and
 * feature edges of the mesh and decimates vertices given
 * by a vertex delete list. */
class MeshDecimation
{
  public:
    typedef std::vector<bool> DeleteList;
    typedef MeshCleanup::VertexIndexRelocList RelocList;

  private:
    typedef std::pair<std::size_t, std::size_t> SplitLine;

  private:
    TriangleMeshPtr mesh;
    FeatureEdgesPtr features;
    VertexInfoListPtr vinfo;
    std::size_t vertex_amount;
    std::size_t vertex_budget;
    DeleteList* dlist;
    RelocList* reloc;

  private:
    bool delete_vertex (std::size_t index);
    void replace_faces (VertexInfo::FaceList aflist,
        MeshFaceList const& flist);

  public:
    MeshDecimation (void);

    /* Sets the mesh and features to be decimated (in-place decimation). */
    void set_mesh (TriangleMeshPtr mesh);
    void set_features (FeatureEdgesPtr features);

    /* Sets the list of vertices to delete. This list may be altered. */
    void set_delete_list (DeleteList& dlist);

    /* Sets an exact vertex budget that the decimation tries to match. */
    void set_exact_budget (std::size_t vertex_amount);

    /* Request a relocation list mapping new to old indices. */
    void fill_vertex_reloc_list (RelocList& reloc);

    /* Run the decimation. */
    void start_decimation (void);

    /* Acquire the final result (note: decimation is in-place). */
    TriangleMeshPtr get_mesh (void) const;
    FeatureEdgesPtr get_features (void) const;
};

/* ---------------------------------------------------------------- */

inline
MeshDecimation::MeshDecimation (void)
{
  this->vertex_budget = 0;
  this->dlist = 0;
  this->reloc = 0;
}

inline void
MeshDecimation::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
MeshDecimation::set_features (FeatureEdgesPtr features)
{
  this->features = features;
}

inline void
MeshDecimation::set_delete_list (MeshDecimation::DeleteList& dlist)
{
  this->dlist = &dlist;
}

inline void
MeshDecimation::set_exact_budget (std::size_t vertex_amount)
{
  this->vertex_budget = vertex_amount;
}

inline void
MeshDecimation::fill_vertex_reloc_list (MeshDecimation::RelocList& reloc)
{
  this->reloc = &reloc;
}


inline TriangleMeshPtr
MeshDecimation::get_mesh (void) const
{
  return this->mesh;
}

inline FeatureEdgesPtr
MeshDecimation::get_features (void) const
{
  return this->features;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_DECIMATION_HEADER */
