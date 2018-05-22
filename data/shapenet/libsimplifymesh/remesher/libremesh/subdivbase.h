#ifndef REMESHER_SUBDIV_BASE_HEADER
#define REMESHER_SUBDIV_BASE_HEADER

#include <vector>

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

class SubdivBase
{
  protected:
    typedef std::vector<std::size_t> NewVertexIndices;

  protected:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    FeatureEdgesPtr features;

  protected:
    virtual void subdiv_impl (void) = 0;
    void insert_faces (NewVertexIndices const& newvi);

  public:
    /* Information for subdivision. All required. */
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    void set_feature_edges (FeatureEdgesPtr features);

    /* Start the algorithm. */
    virtual void start_subdiv (void);
};

/* ---------------------------------------------------------------- */

inline void
SubdivBase::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
SubdivBase::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
SubdivBase::set_feature_edges (FeatureEdgesPtr features)
{
  this->features = features;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_SUBDIV_BASE_HEADER */
