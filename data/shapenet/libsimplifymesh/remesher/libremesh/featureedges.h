#ifndef REMESHER_FEATURE_EDGES_HEADER
#define REMESHER_FEATURE_EDGES_HEADER

#include <vector>

#include "defines.h"
#include "refptrarray.h"
#include "vertexinfo.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

struct FeatureEdgesConf
{
  /* Defines the dihedral angle threashold between faces (in RAD).
   * Any edge with a larger or equal corresponding dihedral angle in
   * the mesh is recognized as feature. (Angle specified in cos(angle)). */
  bool use_angle;
  float angle;

  /* Defines if the mesh boundary should become feature edges. */
  bool border_edges;

  /* Defines if complex edges should become feature edges. */
  bool complex_edges;

  /* Extract contour edges. Just for fun. Angles in RAD. */
  bool extract_contours;
  float angle_theta;
  float angle_phi;

  FeatureEdgesConf (void);
};

/* ---------------------------------------------------------------- */

/*
 * The features edges are extracted using the dilatheral angle of two
 * incident face normals. If this dilatheral angle is larger than the
 * configured threshold, the edge between the faces is marked as feature.
 *
 * FIXME: Some typically irrelevant problems near complex vertices.
 */

class FeatureEdges;
typedef std::vector<std::size_t> FeatureVertexEdges;
typedef RefPtrArray<FeatureEdges, FeatureVertexEdges> FeatureEdgesPtr;

class FeatureEdges : public std::vector<FeatureVertexEdges>
{
  public:
    typedef std::vector<std::size_t> RelocList;

  private:
    FeatureEdgesConf config;
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;

  private:
    bool handle_edge (std::size_t face, std::size_t v1, std::size_t v2);

  protected:
    FeatureEdges (void);

  public:
    static FeatureEdgesPtr create (void);
    FeatureEdgesPtr create_copy (void) const;

    /* Set algorithm configuration. */
    void set_config (FeatureEdgesConf const& config);

    /* Set algorithm data (all required). */
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    TriangleMeshPtr get_mesh (void) const;

    /* Start algorithm; extract all edges matching the config. */
    void extract_features (void);

    /* Feature information. */
    bool is_feature_edge (std::size_t index1, std::size_t index2) const;
    void add_feature_edge (std::size_t index1, std::size_t index2);
    bool rm_feature_edge (std::size_t index1, std::size_t index2);

    FeatureVertexEdges get_edges_for_vertex (std::size_t index) const;
    FeatureVertexEdges get_edges_for_vertex (std::size_t index,
        VertexInfo::VertexList const& vlist) const;
    std::size_t count_edges_for_vertex (std::size_t index,
        VertexInfo::VertexList const& vlist) const;

    /* Fixes feature information if vertex indices have been changed. */
    void fix_index_relocation (FeatureEdges::RelocList const& reloc);
    
    /* Returns the approximate amount of memory used. */
    std::size_t get_memory_usage (void) const;

    /* Dump features to stdout for debugging purposes. */
    void debug_dump (void) const;
};

/* ---------------------------------------------------------------- */

inline
FeatureEdgesConf::FeatureEdgesConf (void)
{
  this->use_angle = true;
  this->angle = std::cos((float)MY_DEG2RAD(40.0f));
  this->border_edges = true;
  this->complex_edges = true;

  this->extract_contours = false;
  this->angle_theta = 0.0f;
  this->angle_phi = 0.0f;
}

inline FeatureEdgesPtr
FeatureEdges::create (void)
{
  return FeatureEdgesPtr(new FeatureEdges);
}

inline
FeatureEdges::FeatureEdges (void)
{
}

inline FeatureEdgesPtr
FeatureEdges::create_copy (void) const
{
  return FeatureEdgesPtr(new FeatureEdges(*this));
}

inline void
FeatureEdges::set_config (FeatureEdgesConf const& config)
{
  this->config = config;
}

inline void
FeatureEdges::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
FeatureEdges::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline TriangleMeshPtr
FeatureEdges::get_mesh (void) const
{
  return this->mesh;
}

inline std::size_t
FeatureEdges::get_memory_usage (void) const
{
  /* This is for the data structure only without features set.
   * Typically, features don't take up much memory. */
  return this->capacity() * sizeof(FeatureVertexEdges);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_FEATURE_EDGES_HEADER */
