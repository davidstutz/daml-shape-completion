#ifndef REMESHER_EDGE_FLIPS_HEADER
#define REMESHER_EDGE_FLIPS_HEADER

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class acts as base class for edge flipping algorithms.
 * It iterates over all edges (but no boundary edges) once and
 * calls a polymorph method to actually check if a flip is to
 * be performed. This is repeated until no more edge flips are
 * performed, or a iteration maximum is reached.
 */
class EdgeFlipsBase
{
  protected:
    std::size_t max_iter;
    std::size_t max_const_flips_iter;
    std::size_t const_flips_iter;
    std::size_t last_flips;

    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    FeatureEdgesPtr features;

  protected:
    /* Iterates over all edges and passes the edge to check_flip_edge
     * for edge flip tests. If the edge is OK to be flipped, the
     * check_flip_edge method is called. */
    std::size_t check_flip_edges (void);

    /* This method must be overwritten. The arguments contain valid
     * vertices and faces, and the edge is safe to be flipped.
     * If the method should return true if the flip should be
     * performed, or return false otherweise, skipping the edge. */
    virtual bool check_flip_edge (std::size_t v1idx,
        std::size_t v2idx, std::size_t v0idx, std::size_t v3idx,
        std::size_t face1, std::size_t face2) = 0;

    /* Actually performs the edge flip updating the data structures. */
    void edge_flip (std::size_t v1idx, std::size_t v2idx,
        std::size_t v0idx, std::size_t v3idx,
        std::size_t face1, std::size_t face2);

  public:
    EdgeFlipsBase (void);

    void set_max_iterations (std::size_t max_iter);
    void set_max_const_flips_iterations (std::size_t iterations);
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    void set_feature_edges (FeatureEdgesPtr features);

    void flip_edges (void);
};

/* ---------------------------------------------------------------- */

inline
EdgeFlipsBase::EdgeFlipsBase (void)
{
  this->max_iter = 0;
  this->max_const_flips_iter = 0;
}

inline void
EdgeFlipsBase::set_max_iterations (std::size_t max_iter)
{
  this->max_iter = max_iter;
}

inline void
EdgeFlipsBase::set_max_const_flips_iterations (std::size_t iterations)
{
  this->max_const_flips_iter = iterations;
}

inline void
EdgeFlipsBase::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
EdgeFlipsBase::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
EdgeFlipsBase::set_feature_edges (FeatureEdgesPtr features)
{
  this->features = features;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_EDGE_FLIPS_HEADER */
