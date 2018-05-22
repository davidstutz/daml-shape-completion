#ifndef REMESHER_RELAXATION_HEADER
#define REMESHER_RELAXATION_HEADER

#include "defines.h"
#include "vec2.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "densityfield.h"
#include "featureedges.h"
#include "vertexref.h"
#include "relocation.h"
#include "microdelaunay.h"
#include "thread.h"

REMESHER_NAMESPACE_BEGIN

struct RelaxationConf
{
  /* Amount of iterations of the relaxation. */
  std::size_t iterations;
  bool pn_lifting;

  RelaxationConf (void);
};

/* ---------------------------------------------------------------- */

class RelaxationThread;
class Relaxation
{
  friend class RelaxationThread;
  protected:
    RelaxationConf config;
    TriangleMeshPtr rmesh;
    TriangleMeshPtr emesh;
    VertexInfoListPtr rvinfo;
    VertexInfoListPtr evinfo;
    DensityFieldPtr rdens;
    FeatureEdgesPtr features;
    FeatureEdgesPtr rfeatures;
    Relocation reloc;

    VertexRefListPtr reflist;
    VertexRefListPtr reflist_new;

  protected:
    /* Process all vertices. */
    void handle_all_vertices (void);
    /* Process all vertices in several threads. */
    void handle_all_vertices_parallel (void);
    /* Processes vertices including start until excluding end. */
    void handle_vertices (std::size_t start, std::size_t end);
    /* Processes a single vertex. */
    void handle_vertex (std::size_t index);

    /* This function queries the new vertex position. */
    virtual Vec2f micropatch_relocate (MicroDelaunay const& md) = 0;

  public:
    Relaxation (void);

    /* Set algorithm configuration. */
    void set_config (RelaxationConf const& config);

    /* Set algorithm data (all required). */
    void set_reference_mesh (TriangleMeshPtr rmesh);
    void set_reference_info (VertexInfoListPtr rvinfo);
    void set_evolving_mesh (TriangleMeshPtr emesh);
    void set_evolving_info (VertexInfoListPtr evinfo);
    void set_evolving_features (FeatureEdgesPtr features);
    void set_reference_features (FeatureEdgesPtr rfeatures);
    void set_vertex_reflist (VertexRefListPtr reflist);
    void set_reference_density (DensityFieldPtr refdens);

    /* Start the algorithm. The base class will do most of the work
     * but querying the new vertex positions on the the new mesh. */
    void start_relaxation (void);
};

/* ---------------------------------------------------------------- */

class RelaxationThread : public Thread
{
  private:
    Relaxation* relax;
    std::size_t start;
    std::size_t end;

  private:
    void* run (void);

  public:
    RelaxationThread (Relaxation* relax, std::size_t start, std::size_t end);
};


/* ---------------------------------------------------------------- */

inline
RelaxationConf::RelaxationConf (void)
{
  this->iterations = 20;
  this->pn_lifting = false;
}

inline
Relaxation::Relaxation (void)
{
}

inline void
Relaxation::set_config (RelaxationConf const& config)
{
  this->config = config;
}

inline void
Relaxation::set_reference_mesh (TriangleMeshPtr rmesh)
{
  this->rmesh = rmesh;
}

inline void
Relaxation::set_evolving_mesh (TriangleMeshPtr emesh)
{
  this->emesh = emesh;
}

inline void
Relaxation::set_reference_info (VertexInfoListPtr rvinfo)
{
  this->rvinfo = rvinfo;
}

inline void
Relaxation::set_evolving_info (VertexInfoListPtr evinfo)
{
  this->evinfo = evinfo;
}

inline void
Relaxation::set_evolving_features (FeatureEdgesPtr features)
{
  this->features = features;
}

inline void
Relaxation::set_reference_features (FeatureEdgesPtr rfeatures)
{
  this->rfeatures = rfeatures;
}

inline void
Relaxation::set_vertex_reflist (VertexRefListPtr reflist)
{
  this->reflist = reflist;
}

inline void
Relaxation::set_reference_density (DensityFieldPtr refdens)
{
  this->rdens = refdens;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_RELAXATION_HEADER */
