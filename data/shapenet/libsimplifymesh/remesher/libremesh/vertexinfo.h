#ifndef REMESHER_VERTEX_INFO_HEADER
#define REMESHER_VERTEX_INFO_HEADER

#include <vector>

#include "defines.h"
#include "refptrarray.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

/*
 * Each vertex is classified into one of SIMPLE, COMPLEX, BORDER
 * and UNREFERENCED. A simple vertex has a full fan of adjacent
 * triangles. A border vertex has a single but incomplete fan
 * of adjacent triangles. An unreferenced vertex has no adjacent
 * triangles. Everything else is a complex vertex, which is
 * basically a non-2-manifold configuration.
 */
enum VertexClass
{
  VERTEX_CLASS_SIMPLE,
  VERTEX_CLASS_COMPLEX,
  VERTEX_CLASS_BORDER,
  VERTEX_CLASS_UNREFERENCED
};

/* ---------------------------------------------------------------- */

/*
 * This class holds per-vertex information.
 *   - The vertex class (see above)
 *   - The adjacent faces (see below)
 */
struct VertexInfo
{
  typedef std::vector<std::size_t> FaceList;
  typedef std::vector<std::size_t> VertexList;

  FaceList adj_faces;
  VertexClass vclass;
};

/* ---------------------------------------------------------------- */

/*
 * This class extracts per-vertex information. Each vertex is cassified
 * into one of the classes SIMPLE, COMPLEX, BORDER and UNREF, see above.
 * Adjacent faces for each vertex are collected and placed into
 * the adj_faces vector. The faces are ordered properly, i.e. in
 * counter-clockwise direction around the vertex. For a boundary vertex,
 * the first and last face in adj_faces are boundary faces.
 * Ordering is not provided for COMPLEX vertices.
 */

class VertexInfoList;
typedef RefPtrArray<VertexInfoList, VertexInfo> VertexInfoListPtr;

class VertexInfoList : public std::vector<VertexInfo>
{
  public:
    typedef std::vector<std::size_t> FaceList;

  private:
    TriangleMeshPtr mesh;

  protected:
    VertexInfoList (void);
    VertexInfoList (TriangleMeshPtr mesh);

  public:
    /* Creates an empty vertex information list. */
    static VertexInfoListPtr create (void);

    /* Create vertex information for the given mesh. */
    static VertexInfoListPtr create (TriangleMeshPtr mesh);

    /* Creates a copy of the vertex information.
     * The new copy will still point to the same mesh. */
    VertexInfoListPtr create_copy (void) const;

    /* Sets the mesh in question and calculates the vertex info. */
    void calc_for_mesh (TriangleMeshPtr mesh);

    /* Sets the mesh in question and clears previous information.
     * The data structure and the classification needs to be
     * manually created if this method is used. */
    void set_mesh (TriangleMeshPtr mesh);

    /* Responsible for preparing the data structure. This
     * need to be called prior any operations on the vertex
     * info and after a mesh has been set. */
    void create_data_structure (void);

    /* Calculates the vertex information for all vertices. */
    void order_and_classify_all (void);

    /* Classifies the vertex referenced by index. As an important
     * byproduct, the adjacent triangles are ordered to form a loop. */
    void order_and_classify (std::size_t index);

    /* Returns a list of adjacent vertices ordered around the vertex. */
    void get_adjacent_vertices (std::size_t index,
        VertexInfo::VertexList& list) const;

    /* Removes an adjacent face from the vertex. Returns true if this
     * face has been removed. Removes one face if it has been found,
     * even if it is present more than once (invalid triangles). */
    bool remove_adjacent_face (std::size_t index, std::size_t face_id);

    /* For an edge given by v1 and v2, the algorithm finds face1 and
     * face2 adjacent to edge v1 -> v2. face1 contains edge v1 -> v2
     * and face 2 contains edge v2 -> v1. v0 and v3 are also found,
     * v0 is the third vertex for face1, v3 is the third vertex for face2.
     * If face1 and face2 is not unique, the method returns false. */
    bool get_edge_info (std::size_t v1, std::size_t v2,
        std::size_t& face1, std::size_t& face2,
        std::size_t& v0, std::size_t& v3);

    /* For a given edge v1->v2, the method returns any face
     * that contains v1->v2 or v2->v1. In the latter case, reverse
     * is set to true. The edge offset within the face is also returned.
     * If v1->v2 is not a valid edge, face and offset is set to MAX_SIZE_T.
     * Any passed pointer may be NULL to omit that information. */
    void get_face_for_edge (std::size_t v1, std::size_t v2,
        std::size_t* face, std::size_t* edge_off, bool* reverse);

    /* Returns true iff v1 -> v2 is an edge in the mesh. */
    bool is_mesh_edge (std::size_t v1, std::size_t v2) const;

    /* Returns the list of faces that share the given edge. */
    FaceList get_faces_for_edge (std::size_t v1, std::size_t v2);

    /* Checks if the given loop of vertices is complex. A loop is
     * complex if two vertices that are non-neighboring in the loop
     * are connected with an edge in the mesh. Thus the two vertices
     * ARE neighbors but not in the loop. Complex loops should not be
     * simplified, otherwise topological changes may be introduced. */
    bool is_complex_vertex_loop (VertexInfo::VertexList const& vlist);

    /* Returns the amount of memory for the data structure. */
    std::size_t get_memory_usage (void) const;

    void debug_vertex (std::size_t vid) const;
};

/* ---------------------------------------------------------------- */

inline
VertexInfoList::VertexInfoList (void)
{
}

inline
VertexInfoList::VertexInfoList (TriangleMeshPtr mesh)
{
  this->calc_for_mesh(mesh);
}

inline VertexInfoListPtr
VertexInfoList::create (void)
{
  return VertexInfoListPtr(new VertexInfoList);
}

inline VertexInfoListPtr
VertexInfoList::create (TriangleMeshPtr mesh)
{
  return VertexInfoListPtr(new VertexInfoList(mesh));
}

inline VertexInfoListPtr
VertexInfoList::create_copy (void) const
{
  return VertexInfoListPtr(new VertexInfoList(*this));
}

inline void
VertexInfoList::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
VertexInfoList::calc_for_mesh (TriangleMeshPtr mesh)
{
  this->set_mesh(mesh);
  this->create_data_structure();
  this->order_and_classify_all();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_VERTEX_INFO_HEADER */
