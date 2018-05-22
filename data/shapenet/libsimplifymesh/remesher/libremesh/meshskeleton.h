#ifndef REMESHER_MESH_SKELETON_HEADER
#define REMESHER_MESH_SKELETON_HEADER

#include <list>
#include <vector>

#include "defines.h"
#include "refptr.h"
#include "featureedges.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This structure represents a backbone. It may be either closed
 * or open, and it contains all vertices of the backbone in order.
 */
struct MeshBackbone
{
  typedef std::list<std::size_t> Vertices;
  typedef Vertices::const_iterator Iter;

  bool closed;
  Vertices verts;
};

typedef std::list<std::size_t>::const_iterator BackboneIter;

/* ---------------------------------------------------------------- */

/*
 * This class extracts the feature skeleton of a mesh for sequential access.
 * It basically takes the default feature represenation and converts
 * it into data structures better suited for resampling.
 *
 * Two different types of backbones are recognized: Closed backbones which
 * are cyclic, and open backbones between two corner vertices.
 */

class MeshSkeleton;
typedef RefPtrArray<MeshSkeleton, MeshBackbone> MeshSkeletonPtr;

class MeshSkeleton : public std::vector<MeshBackbone>
{
  protected:
    typedef std::pair<std::size_t, std::size_t> FeatureEdge;

  private:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    FeatureEdgesPtr features;
    float max_angle;
    std::size_t num_closed;
    std::size_t num_open;

  protected:
    MeshSkeleton (void);

    /* Finds the next edge for the given feature edge. The function returns
     * true on success, or false if advancing is not possible. */
    bool advance_edge (FeatureEdge& edge) const;

    /* Starts collecting from a given seed edge.
     * Removes consumed edges from the pool data structure. */
    void collect (FeatureEdge edge, FeatureEdgesPtr pool);

  public:
    static MeshSkeletonPtr create (void);

    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    void set_feature_edges (FeatureEdgesPtr features);

    /* Sets a maximum angle between two subsequent edges on a feature
     * crease. If the angle is exceeded, a corner is inserted. */
    void set_max_angle (float angle);

    void extract (void);

    /* Returns the amount of closed backbones. */
    std::size_t get_closed_amount (void) const;
    /* Returns the amount of open backbones. */
    std::size_t get_open_amount (void) const;
    /* Returns the amount of corner vertices. */
    std::size_t get_corner_amount (void) const;

    void debug_print (void) const;
};

/* ---------------------------------------------------------------- */

inline
MeshSkeleton::MeshSkeleton (void)
{
  this->num_closed = 0;
  this->num_open = 0;
  this->max_angle = 0.0f;
}

inline MeshSkeletonPtr
MeshSkeleton::create (void)
{
  return MeshSkeletonPtr(new MeshSkeleton);
}

inline void
MeshSkeleton::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
MeshSkeleton::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
MeshSkeleton::set_feature_edges (FeatureEdgesPtr features)
{
  this->features = features;
}

inline void
MeshSkeleton::set_max_angle (float angle)
{
  this->max_angle = angle;
}

inline std::size_t
MeshSkeleton::get_closed_amount (void) const
{
  return this->num_closed;
}

inline std::size_t
MeshSkeleton::get_open_amount (void) const
{
  return this->num_open;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_SKELETON_HEADER */
