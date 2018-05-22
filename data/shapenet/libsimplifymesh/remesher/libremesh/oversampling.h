#ifndef REMESHER_OVERSAMPLING_HEADER
#define REMESHER_OVERSAMPLING_HEADER

#include "defines.h"
#include "trianglemesh.h"
#include "densityfield.h"
#include "vertexinfo.h"
#include "vertexref.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

struct OversamplingConf
{
  bool use_loop_subdiv;
  bool use_linear_subdiv;

  bool use_slicing;
  bool lazy_triangulation;
  std::size_t slices_x;
  std::size_t slices_y;
  std::size_t slices_z;

  OversamplingConf (void);
};

/* ---------------------------------------------------------------- */

class Oversampling
{
  private:
    OversamplingConf config;

    TriangleMeshPtr mesh;
    DensityFieldPtr density;
    FeatureEdgesPtr features;
    VertexInfoListPtr vinfo;

  protected:
    void run_loop_subdiv (void);
    void run_linear_subdiv (void);
    void run_grid_slicing (void);

  public:
    Oversampling (void);

    void set_config (OversamplingConf const& conf);
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    void set_density_field (DensityFieldPtr density);
    void set_feature_edges (FeatureEdgesPtr features);

    void start_oversampling (void);
};

/* ---------------------------------------------------------------- */

inline
OversamplingConf::OversamplingConf (void)
{
  this->use_loop_subdiv = true;
  this->use_linear_subdiv = false;

  this->use_slicing = false;
  this->lazy_triangulation = true;
  this->slices_x = 10;
  this->slices_y = 10;
  this->slices_z = 10;
}

inline
Oversampling::Oversampling (void)
{
}

inline void
Oversampling::set_config (OversamplingConf const& conf)
{
  this->config = conf;
}

inline void
Oversampling::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
Oversampling::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
Oversampling::set_density_field (DensityFieldPtr density)
{
  this->density = density;
}

inline void
Oversampling::set_feature_edges (FeatureEdgesPtr features)
{
  this->features = features;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_OVERSAMPLING_HEADER */
