#ifndef REMESHER_TRIANGULATOR_HEADER
#define REMESHER_TRIANGULATOR_HEADER

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class is capable of triangulating a set of vertices given
 * as vertex loop, i.e. a polygon of vertices. The vertices are
 * indices to a reference mesh. The 3D positions of the vertices
 * to calculate the triangulation is taken from the mesh.
 *
 * The method in question is a recursive loop-splitting procedure
 * that is used in a simplification algorithm described in:
 *
 *   William J. Schroeder, Jonathan A. Zarge, William E. Lorensen
 *   Decimation of Triangle Meshes
 *
 * An initial splitline can be specified, if there is a constrain
 * on a certain edge (i.e. a feature edge that is simplified).
 */
class Triangulator
{
  public:
    typedef std::pair<std::size_t, std::size_t> SplitLine;

  protected:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;

  protected:
    /* Checks the quality of the split line. If the split line is
     * invalid, quality is negative. Otherwise quality is positive. */
    float check_splitline_quality (SplitLine const& sl,
        VertexInfo::VertexList const& vl, Vec3f const& ap_normal);

  public:
    Triangulator (void);
    Triangulator (TriangleMeshPtr mesh);

    void set_mesh (TriangleMeshPtr mesh);

    /* Setting the vertex info enables checks to avoid "invalid"
     * triangulations (those that change genus and so on) to happen. */
    void set_vertex_info (VertexInfoListPtr vinfo);

    /* Meshes the given loop of vertices by recursive loop splitting.
     * The overlapping check needs the normal of the average plane,
     * the normal may not be normalized. A splitline may be specified
     * for the first split. */
    bool triangulate (VertexInfo::VertexList const& vlist,
        MeshFaceList& flist, Vec3f const& ap_normal,
        SplitLine const& sline);
};

/* ---------------------------------------------------------------- */

inline
Triangulator::Triangulator (void)
{
}

inline
Triangulator::Triangulator (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
Triangulator::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
Triangulator::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_TRIANGULATOR_HEADER */
