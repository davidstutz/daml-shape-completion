#ifndef REMESHER_DELAUNAY_FLIPS_HEADER
#define REMESHER_DELAUNAY_FLIPS_HEADER

#include "defines.h"
#include "edgeflips.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class performs edge flips directly on the 3D mesh to
 * obtain a delaunay triangulation of the mesh.
 * Edge flips are performed if the sum of the angles opposite
 * to the edge sum up to less than 180 degrees.
 */
class DelaunayFlips : public EdgeFlipsBase
{
  protected:
    bool check_flip_edge (std::size_t v1idx,
        std::size_t v2idx, std::size_t v0idx, std::size_t v3idx,
        std::size_t face1, std::size_t face2);

  public:
    DelaunayFlips (void);
    DelaunayFlips (TriangleMeshPtr mesh, VertexInfoListPtr vinfo);
};

/* ---------------------------------------------------------------- */

inline
DelaunayFlips::DelaunayFlips (void)
{
  this->set_max_iterations(30);
  this->set_max_const_flips_iterations(3);
}

inline
DelaunayFlips::DelaunayFlips (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
{
  this->set_mesh(mesh);
  this->set_vertex_info(vinfo);
  this->set_max_iterations(30);
  this->set_max_const_flips_iterations(3);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_DELAUNAY_FLIPS_HEADER */
