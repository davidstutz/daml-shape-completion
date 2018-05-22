#ifndef REMESHER_INTERFACE_HEADER
#define REMESHER_INTERFACE_HEADER

#include <string>

#include "defines.h"
#include "vertexinfo.h"
#include "vertexref.h"
#include "trianglemesh.h"
#include "simplification.h"
#include "densityfield.h"
#include "featureedges.h"
#include "oversampling.h"
#include "resampling.h"
#include "relaxlloyd.h"
#include "relaxarea.h"

REMESHER_NAMESPACE_BEGIN

class Interface
{
  private:
    /* Reference mesh data. */
    TriangleMeshPtr reference_mesh;
    VertexInfoListPtr reference_vinfo;
    DensityFieldPtr reference_density;
    FeatureEdgesPtr reference_features;

    /* Evolving mesh data. */
    TriangleMeshPtr evolving_mesh;
    VertexInfoListPtr evolving_vinfo;
    VertexRefListPtr evolving_vref;
    FeatureEdgesPtr evolving_features;

    /* Algorithm configurations. */
    OversamplingConf oversampling_conf;
    DensityFieldConf density_field_conf;
    FeatureEdgesConf feature_edges_conf;
    SimplificationConf simplification_conf;
    ResamplingConf resampling_conf;
    RelaxationConf lloyd_conf;
    RelaxationConf area_equal_conf;

  protected:
    void assure_reference_mesh (void);
    void assure_reference_vinfo (void);
    void assure_reference_features (void);
    void assure_evolving_mesh (void);
    void assure_vertex_references (void);

  public:
    Interface (void);

    /* Provides access to the models. */
    TriangleMeshPtr get_reference_mesh (void) const;
    TriangleMeshPtr get_evolving_mesh (void) const;

    FeatureEdgesPtr get_reference_features (void);
    FeatureEdgesPtr get_evolving_features (void) const;

    /* Loads the reference mesh. */
    void load_model (std::string const& modelfile);

    /* Configuring and executing oversampling. */
    void set_oversampling_conf (OversamplingConf const& conf);
    void exec_oversampling (void);

    /* Configuring and executing the feature extracation. */
    void set_feature_edges_conf (FeatureEdgesConf const& conf);
    void exec_feature_extraction (void);
    void clear_feature_edges (void);

    /* Configuring and executing the density calculation. */
    void set_density_field_conf (DensityFieldConf const& conf);
    void exec_density_calculation (void);
    void clear_density_field (void);
    DensityFieldPtr get_density_field (void) const;

    /* Configuring and executing the simplification. */
    void set_simplification_conf (SimplificationConf const& conf);
    void exec_simplification (void);

    /* Execute resampling of the mesh. */
    void set_resampling_conf (ResamplingConf const& conf);
    void exec_resampling (void);

    /* Executing the mesh cleanup. */
    void clean_reference_mesh (void);

    /* Remove duplicated vertices. */
    void clean_duplicated_vertices (void);

    /* Copies the reference mesh to the evolving mesh. */
    void copy_reference_mesh (void);

    /* Copies the evolving mesh to the reference mesh. */
    void copy_evolving_mesh (void);

    /* Some information about the evolving mesh. */
    VertexRefListPtr get_evolving_vertex_refs (void) const;
    VertexInfoListPtr get_evolving_vertex_info (void) const;

    /* Optimizes evolving/reference mesh and adapts all evolving data. */
    void optimize_reference_mesh (void);
    void optimize_evolving_mesh (void);

    /* Configuring and executing the Lloyd relaxation. */
    void set_lloyd_conf (RelaxationConf const& conf);
    void exec_lloyd (void);

    /* Configuring and executing the area equalization. */
    void set_area_equal_conf (RelaxationConf const& conf);
    void exec_area_equalization (void);

    /* Restore Delaunay criterion. */
    void exec_delaunay_flips (void);
};

/* ---------------------------------------------------------------- */

inline TriangleMeshPtr
Interface::get_reference_mesh (void) const
{
  return this->reference_mesh;
}

inline TriangleMeshPtr
Interface::get_evolving_mesh (void) const
{
  return this->evolving_mesh;
}

inline void
Interface::set_density_field_conf (DensityFieldConf const& conf)
{
  this->density_field_conf = conf;
}

inline void
Interface::clear_density_field (void)
{
  this->reference_density->clear();
}

inline DensityFieldPtr
Interface::get_density_field (void) const
{
  return this->reference_density;
}

inline void
Interface::set_feature_edges_conf (FeatureEdgesConf const& conf)
{
  this->feature_edges_conf = conf;
}

inline void
Interface::clear_feature_edges (void)
{
  this->reference_features->clear();
}

inline FeatureEdgesPtr
Interface::get_reference_features (void)
{
  this->assure_reference_features();
  return this->reference_features;
}

inline FeatureEdgesPtr
Interface::get_evolving_features (void) const
{
  return this->evolving_features;
}

inline VertexRefListPtr
Interface::get_evolving_vertex_refs (void) const
{
  return this->evolving_vref;
}

inline VertexInfoListPtr
Interface::get_evolving_vertex_info (void) const
{
  return this->evolving_vinfo;
}

inline void
Interface::set_simplification_conf (SimplificationConf const& conf)
{
  this->simplification_conf = conf;
}

inline void
Interface::set_resampling_conf (ResamplingConf const& conf)
{
  this->resampling_conf = conf;
}

inline void
Interface::set_oversampling_conf (OversamplingConf const& conf)
{
  this->oversampling_conf = conf;
}

inline void
Interface::set_lloyd_conf (RelaxationConf const& conf)
{
  this->lloyd_conf = conf;
}

inline void
Interface::set_area_equal_conf (RelaxationConf const& conf)
{
  this->area_equal_conf = conf;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_INTERFACE_HEADER */
