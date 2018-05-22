#ifndef REMESHER_MESH_STATS_HEADER
#define REMESHER_MESH_STATS_HEADER

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

class MeshStats
{
  private:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;

  public:
    MeshStats (TriangleMeshPtr mesh, VertexInfoListPtr vinfo);
    void print_stats (void);
};

/* ---------------------------------------------------------------- */

inline
MeshStats::MeshStats (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
  : mesh(mesh), vinfo(vinfo)
{
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_STATS_HEADER */
