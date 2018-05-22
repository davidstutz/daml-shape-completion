#ifndef REMESHER_PATCH_3D_HEADER
#define REMESHER_PATCH_3D_HEADER

#include <set>
#include <string>

#include "refptr.h"
#include "vertexinfo.h"
#include "trianglemesh.h"
#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/*
 * The classes in this file are capable of creating triangle patches
 * in three dimensions. Different strategies on how to collect the
 * triangles can be implemented. If the patch is to be flattened
 * to a circular boundary, the 3d patch should similarly be circular
 * for optimal results.
 */

class Patch3d;
typedef RefPtr<Patch3d> Patch3dPtr;

class Patch3d
{
  public:
    typedef std::set<std::size_t> MeshFaceSet;

  protected:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    MeshFaceSet faces;

  protected:
    Patch3d (void);
    void set_data (TriangleMeshPtr mesh, VertexInfoListPtr vinfo);

  public:
    TriangleMeshPtr get_mesh (void) const;
    MeshFaceSet const& get_face_set (void) const;

    /* This is inefficient! Stores all vertices and only a subset of faces. */
    void write_debug_mesh (std::string const& filename);
};

/* ---------------------------------------------------------------- */

/*
 * This class is capable of creating a dynamic patch of triangles
 * from a 3D mesh. The patch is created using a method from
 *
 *   Vitaly Surazhsky and Craig Gotsman
 *   Explicit Surface Remeshing
 *
 * A patch is specified with three face indices on a reference mesh.
 * Triangles are collected in a breadth-first search over the faces
 * until the three given triangles are visited. The resulting
 * patch is circular and contails the "triangle" specified by the
 * three faces.
 */

class Patch3dTriangle;
typedef RefPtr<Patch3dTriangle> Patch3dTrianglePtr;

class Patch3dTriangle : public Patch3d
{
  protected:
    Patch3dTriangle (void);

    void add_faces (std::size_t f1, std::size_t f2, std::size_t f3);
    bool add_face (std::size_t f);

  public:
    static Patch3dTrianglePtr create (TriangleMeshPtr mesh,
        VertexInfoListPtr vinfo);

    void create_patch (std::size_t f1, std::size_t f2, std::size_t f3);
};

/* ---------------------------------------------------------------- */

/*
 * Simple patch class that creates a 3D patch for a vertex and
 * a certain neighborhood (1-ring, 2-ring, etc) of the vertex.
 */
class Patch3dNRing;
typedef RefPtr<Patch3dNRing> Patch3dNRingPtr;

class Patch3dNRing : public Patch3d
{
  public:
    static Patch3dNRingPtr create (TriangleMeshPtr mesh,
        VertexInfoListPtr vinfo);

    void create_patch (std::size_t index, std::size_t n_rings);
};

/* ---------------------------------------------------------------- */

inline
Patch3d::Patch3d (void)
{
}

inline void
Patch3d::set_data (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
{
  this->mesh = mesh;
  this->vinfo = vinfo;
}

inline TriangleMeshPtr
Patch3d::get_mesh (void) const
{
  return this->mesh;
}

inline Patch3d::MeshFaceSet const&
Patch3d::get_face_set (void) const
{
  return this->faces;
}

inline
Patch3dTriangle::Patch3dTriangle (void)
{
}

inline Patch3dTrianglePtr
Patch3dTriangle::create (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
{
  Patch3dTrianglePtr ret(new Patch3dTriangle);
  ret->set_data(mesh, vinfo);
  return ret;
}

inline Patch3dNRingPtr
Patch3dNRing::create (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
{
  Patch3dNRingPtr ret(new Patch3dNRing);
  ret->set_data(mesh, vinfo);
  return ret;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PATCH_3D_HEADER */
