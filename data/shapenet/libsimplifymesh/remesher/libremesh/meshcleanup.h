#ifndef REMESHER_MESH_CLEANUP_HEADER
#define REMESHER_MESH_CLEANUP_HEADER

#include <vector>

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class cleans up a triangle mesh. Cleanup requires vertex
 * inforation. This vertex information is invalidated during cleanup.
 *
 * It can be used to remove vertices and faces which are marked as deleted.
 * In this case the "cleanup_deleted" method is to be used.
 * This method removes all vertices from the data structure that are marked
 * as deleted in the passed argument. It also removes all triangles
 * indexing vertex 0, i.e. set all three vertex indices to 0 to remove a
 * face. The list of deleted vertices is invalidated during cleanup.
 * Since vector resizing is used after cleanup, the mesh must be copied
 * to actually free up memory as a result from the cleanup.
 *
 * The class can also be used to remove duplicated vertices (that is
 * vertices that are very close to each other, depending on epsilon). The
 * triangulation is fixed and some faces are deleted (common faces for the
 * vertices to be merged). Unreferenced vertices are also removed.
 * Method "cleanup_mesh" is to be used for this feature.
 * Note that close vertices are ONLY removed IF there is a face
 * connecting the two vertices in question.
 *
 * The class offers to supply a vertex index relocation list. After
 * cleanup, this list maps the index of a vertex to its old position
 * and vertices can still be associated with the original mesh.
 * No such strucure is created by default.
 */

class MeshCleanup
{
  public:
    typedef std::vector<std::size_t> VertexIndexRelocList;
    typedef std::vector<bool> VertexDeletedList;

  private:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    VertexIndexRelocList* ireloc;

  private:
    /* Merges two vertices and deletes all zero-area faces resulting
     * from the merge. The amount of removed triangles is returned.
     * If no triangle was removed, the vertices could not be merged,
     * otherwise the vertex `b' is unused should be deleted. */
    std::size_t merge_vertices (std::size_t a, std::size_t b);

    /* Deletes a face and removes this face from adjacent vertices. */
    bool delete_invalid_face (std::size_t face);

  public:
    MeshCleanup (TriangleMeshPtr mesh);

    /* Moving vertices requires fixing faces. The vertex info
     * is neccessary for fast operation of the algorithm. */
    void set_vertex_info (VertexInfoListPtr vinfo);

    /* The vertex relocation list is filled during cleanup.
     * It maps vertex indices to their old positions in the mesh. */
    void set_vertex_reloc_list (VertexIndexRelocList* irlist);
    void unset_vertex_reloc_list (void);

    /*
     * The main algorithms.
     */

    /* Cleanup of deleted vertices. Note that also normals are deleted. */
    void cleanup_deleted (VertexDeletedList& dlist);

    /* Cleans up duplicated (connected) and unreferenced vertices. */
    void cleanup_mesh (float epsilon = MY_FLT_EPS);

    /* Cleans up duplicated (unconnected) vertices. A bit expensive. */
    void cleanup_duplicated (float epsilon = MY_FLT_EPS);
};

/* ---------------------------------------------------------------- */

inline
MeshCleanup::MeshCleanup (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
  this->ireloc = 0;
}

inline void
MeshCleanup::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
MeshCleanup::set_vertex_reloc_list (VertexIndexRelocList* irlist)
{
  this->ireloc = irlist;
}

inline void
MeshCleanup::unset_vertex_reloc_list (void)
{
  this->ireloc = 0;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_CLEANUP_HEADER */
