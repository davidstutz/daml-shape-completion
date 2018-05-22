#ifndef REMESHER_TRIANGLE_MESH_HEADER
#define REMESHER_TRIANGLE_MESH_HEADER

#include <vector>

#include "defines.h"
#include "refptr.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

/*
 * Determines if angle-weighted pseudo normals should be used for
 * vertex normal calculation. Otherwise, normals are calculated
 * by averaging area-weighted adjacent face normals.
 */
#define MESH_AWPN_NORMALS 1

class TriangleMesh;
typedef RefPtr<TriangleMesh> TriangleMeshPtr;

/*
 * Some types that should be used outside this class.
 * Unfortunately, the vertex indices need to be of type "unsigned int"
 * because OpenGL requires it for fast vertex array rendering.
 * Since this might change, it's also typedefed.
 */
typedef unsigned int MeshVIndex;
typedef std::vector<Vec3f> MeshVertexList;
typedef std::vector<MeshVIndex> MeshFaceList;
typedef std::vector<Vec3f> MeshNormalList;
typedef std::vector<Vec3uc> MeshVertexColorList;

/*
 * This class holds a triangle mesh as a list of vertices
 * and a list of faces. The list of faces contains three
 * elements per face indexing the three vertices of the face.
 *
 * The class also cares about creating normals for faces and vertices
 * if requested. Faces can be inverted by inverting the vertex order
 * of each face. For displaying purposes, the model can be scaled to
 * a bounding cube of unit length and centered in the coordinate system.
 */

class TriangleMesh
{
  protected:
    MeshVertexList vertices;
    MeshFaceList faces;
    MeshNormalList face_normals;
    MeshNormalList vertex_normals;
    MeshVertexColorList vertex_colors;

  protected:
    TriangleMesh (void);

  public:
    static TriangleMeshPtr create (void);
    TriangleMeshPtr create_copy (bool with_normals = true) const;

    /* Information about vertices, faces and normals. */
    MeshVertexList const& get_vertices (void) const;
    MeshFaceList const& get_faces (void) const;
    MeshNormalList const& get_face_normals (void) const;
    MeshNormalList const& get_vertex_normals (void) const;
    MeshVertexColorList const& get_vertex_colors (void) const;

    MeshVertexList& get_vertices (void);
    MeshFaceList& get_faces (void);
    MeshNormalList& get_face_normals (void);
    MeshNormalList& get_vertex_normals (void);
    MeshVertexColorList& get_vertex_colors (void);

    /* Calculation of normals. */
    void ensure_normals (bool face = true, bool vertex = true);
    void recalc_normals (bool face = true, bool vertex = true);

    /* Invert orientation of faces. */
    void invert_faces (void);

    /* Scale and center the model in the coordinate system. */
    void scale_and_center (void);

    /* Clear model data. */
    void clear (void);
    void clear_normals (void);

    /* Memory information. */
    void memory_debug (void) const;
    std::size_t get_memory_usage (void) const;
};

/* ---------------------------------------------------------------- */

inline
TriangleMesh::TriangleMesh (void)
{
}

inline TriangleMeshPtr
TriangleMesh::create (void)
{
  return TriangleMeshPtr(new TriangleMesh);
}

inline TriangleMeshPtr
TriangleMesh::create_copy (bool with_normals) const
{
  TriangleMeshPtr ret(new TriangleMesh);
  ret->get_vertices() = this->vertices;
  ret->get_faces() = this->faces;
  if (with_normals)
  {
    ret->get_vertex_normals() = this->vertex_normals;
    ret->get_face_normals() = this->face_normals;
  }
  return ret;
}

inline MeshVertexList const&
TriangleMesh::get_vertices (void) const
{
  return this->vertices;
}

inline MeshFaceList const&
TriangleMesh::get_faces (void) const
{
  return this->faces;
}

inline MeshNormalList const&
TriangleMesh::get_face_normals (void) const
{
  return this->face_normals;
}

inline MeshNormalList const&
TriangleMesh::get_vertex_normals (void) const
{
  return this->vertex_normals;
}

inline MeshVertexColorList const&
TriangleMesh::get_vertex_colors (void) const
{
  return this->vertex_colors;
}

inline MeshVertexList&
TriangleMesh::get_vertices (void)
{
  return this->vertices;
}

inline MeshFaceList&
TriangleMesh::get_faces (void)
{
  return this->faces;
}

inline MeshNormalList&
TriangleMesh::get_face_normals (void)
{
  return this->face_normals;
}

inline MeshNormalList&
TriangleMesh::get_vertex_normals (void)
{
  return this->vertex_normals;
}

inline MeshVertexColorList&
TriangleMesh::get_vertex_colors (void)
{
  return this->vertex_colors;
}

inline void
TriangleMesh::clear_normals (void)
{
  this->vertex_normals.clear();
  this->face_normals.clear();
}

inline void
TriangleMesh::clear (void)
{
  this->vertices.clear();
  this->faces.clear();
  this->face_normals.clear();
  this->vertex_normals.clear();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_TRIANGLE_MESH_HEADER */
