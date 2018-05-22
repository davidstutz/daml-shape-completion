#ifndef REMESHER_RESAMPLING_HEADER
#define REMESHER_RESAMPLING_HEADER

#include <set>
#include <vector>

#include "defines.h"
#include "vertexref.h"
#include "vertexinfo.h"
#include "densityfield.h"
#include "featureedges.h"
#include "trianglemesh.h"
#include "meshskeleton.h"
#include "meshdecimation.h"

#define RESAMPLING_ON_MESHCOPY 1
#define MUTUAL_TESSELLATION_COLORS 1
#define MUTUAL_COLOR_ORIG Vec3uc(255, 0, 0)
#define MUTUAL_COLOR_NEW Vec3uc(0, 0, 255)

REMESHER_NAMESPACE_BEGIN

/*
 * Mesh resampling configuration.
 */
struct ResamplingConf
{
  /* The amount of samples to distribute on the surface and features. */
  std::size_t sample_amount;

  /* When inserting in the original mesh, decimate vertices afterwards. */
  bool perform_decimation;

  ResamplingConf (void);
};

/* ---------------------------------------------------------------- */

/* Simple data structure that groups mesh input data. */
struct InputMeshData
{
  TriangleMeshPtr mesh;
  VertexInfoListPtr vinfo;
  FeatureEdgesPtr features;
  DensityFieldPtr density;
};

/* Simple data structure that groups mesh output data. */
struct OutputMeshData
{
  TriangleMeshPtr mesh;
  FeatureEdgesPtr features;
  VertexRefListPtr vrefs;
};

/* ---------------------------------------------------------------- */

/*
 * A small helper struct that represents an edge (v1->v2) and the
 * amount of samples for that edge. The actual sample positions on
 * the edge are stored elsewhere. Edges with v1 = v2 are corner vertices.
 */
struct SamplesPerEdge
{
  MeshVIndex v1;
  MeshVIndex v2;
  std::size_t samples;

  SamplesPerEdge (MeshVIndex v1, MeshVIndex v2, std::size_t samples)
    : v1(v1), v2(v2), samples(samples) {}
#if __WORDSIZE==64
  SamplesPerEdge (std::size_t v1, std::size_t v2, std::size_t samples)
    : v1((MeshVIndex)v1), v2((MeshVIndex)v2), samples(samples) {}
#endif
};

/* ---------------------------------------------------------------- */

/* Result type for samples per face. */
typedef std::vector<std::size_t> SamplesPerFaceList;

/* Result types for samples on features. */
typedef std::vector<SamplesPerEdge> SamplesPerEdgeList;
typedef std::vector<float> EdgeSamplePositions;

/*
 * The sampler for the resampling procedure. The sampler is first
 * calibrated by analyzing the mesh surface and features. It then
 * distributes the given amount of samples on the surface and features
 * accordingly. The result is given in surface samples per face of the
 * input mesh, and exact sample positions on feature creases.
 */
class ResamplingSampler
{
  private:
    /* Input mesh and data structures. */
    InputMeshData rmesh;

    /* Working data. */
    MeshSkeletonPtr skel;

    /* Output data structures: Samples per face and on features. */
    SamplesPerFaceList per_face;
    SamplesPerEdgeList per_edge;
    EdgeSamplePositions per_edge_pos;

  protected:
    /* Iterates over all faces and sums the face density. */
    void calc_surface_density (float& area, float& density) const;
    //float get_face_density (std::size_t fid) const;

    /* Iterates over all feature eges and sums edge density. */
    void calc_feature_density (float& length, float& density) const;
    float get_backbone_density (MeshBackbone const& mb) const;
    //float get_edge_density (std::size_t vid1, std::size_t vid2) const;

    /* Resamples the corner vertices. */
    std::size_t sample_corners (void);
    /* Resamples all feature creases. */
    std::size_t sample_features (std::size_t samples, float total_density);
    /* Resamples the (smooth) model surface. */
    std::size_t sample_surface (std::size_t samples, float total_density);

    /* Samples the given backbone with the given amount of samples. */
    void sample_backbone (MeshBackbone const& mb, float backbone_dens,
        std::size_t samples);

    /* Samples a single edge. An excess teleport may be specified,
     * this should be positive or zero and causes samples to be placed
     * earlier (after excess_teleport mass has been collected).
     * A new excess teleport is set after sampling.
     * The number of sampled vertices is returned. */
    std::size_t sample_edge (std::size_t v1, std::size_t v2,
        float density_spacing, float& excess_teleport,
        std::size_t max_vertices);

    /* Evaluate new position on an edge for a given density.
     * Returns a new coordinate [0,len] as new position.
     * The edge is specified with a length and density values d1 and d2.
     * The progress to be made is given in density from a given position.
     */
    float eval_edge_coord(float len, float d1, float d2, float cur, float dens);

    /* Place samples in the result data structure. */
    void place_feature_sample (std::size_t v1, std::size_t v2, float lambda);

  public:
    ResamplingSampler (void);

    /* Set the reference mesh and additional data structures. */
    void set_mesh_data (InputMeshData const& rmesh);

    /* Execute the algorithm. */
    void exec_resampling (std::size_t sample_amount);

    /* Sampler output. */
    SamplesPerFaceList const& get_samples_per_face (void) const;
    SamplesPerEdgeList const& get_samples_per_edge (void) const;
    EdgeSamplePositions const& get_edge_sample_pos (void) const;
};

/* ---------------------------------------------------------------- */

/*
 * The mesher takes as input the amount of samples per face and
 * the exact sample positions on feature edges and meshes the result
 * directly in 3D space. This results in a new, resampled mesh.
 */
class ResamplingMesher
{
  private:
    InputMeshData rmesh;
    OutputMeshData emesh;
    ResamplingConf config;

    ResamplingSampler const* sampler;
    MeshDecimation::DeleteList dlist;

  private:
    void sample_faces (SamplesPerFaceList const& flist);
    void sample_features (SamplesPerEdgeList const& elist,
        EdgeSamplePositions const& epos);

    /* Samples the given triangle with the given amount of samples. */
    void sample_triangle (std::size_t tri, std::size_t samples);

    /* Returns a random barycentric coordinate. */
    Vec3f get_random_bary (void) const;

    /* Returns random barycentric coordinate for density distribution. */
    Vec3f get_random_bary (float d1, float d2, float d3);

  public:
    ResamplingMesher (void);

    void set_input_mesh (InputMeshData const& rmesh);
    void set_config (ResamplingConf const& conf);
    void set_sampler (ResamplingSampler const& sampler);

    void exec_meshing (void);

    /* Returns the final mesh. */
    OutputMeshData get_output_mesh (void) const;
};

/* ---------------------------------------------------------------- */

/*
 * Resampling class is a helper class that ties together the sample
 * component and the meshing component to produce the final mesh.
 */
class Resampling
{
  private:
    ResamplingConf config;

    ResamplingSampler sampler;
    ResamplingMesher mesher;

    InputMeshData rmesh;
    OutputMeshData emesh;    

  public:
    Resampling (void);

    /* Set algorithm config. */
    void set_config (ResamplingConf const& config);

    /* Set the reference mesh and additional data structures. */
    void set_reference_mesh (TriangleMeshPtr rmesh);
    void set_reference_vinfo (VertexInfoListPtr rvinfo);
    void set_reference_features (FeatureEdgesPtr rfeatures);
    void set_reference_density (DensityFieldPtr rdensity);

    /* Start the resampling pipeline. */
    void exec_resampling (void);

    /* Result after resampling. */
    TriangleMeshPtr get_resampled_mesh (void) const;
    FeatureEdgesPtr get_resampled_features (void) const;
    VertexRefListPtr get_vertex_references (void) const;
};

/* -------------------- Implementation: Config -------------------- */

inline
ResamplingConf::ResamplingConf (void)
{
  this->sample_amount = 10000;
  this->perform_decimation = true;
}

/* ------------------- Implementation: Sampler -------------------- */

inline
ResamplingSampler::ResamplingSampler (void)
{
}

inline void
ResamplingSampler::set_mesh_data (InputMeshData const& rmesh)
{
  this->rmesh = rmesh;
}

inline SamplesPerFaceList const&
ResamplingSampler::get_samples_per_face (void) const
{
  return this->per_face;
}

inline SamplesPerEdgeList const&
ResamplingSampler::get_samples_per_edge (void) const
{
  return this->per_edge;
}

inline EdgeSamplePositions const&
ResamplingSampler::get_edge_sample_pos (void) const
{
  return this->per_edge_pos;
}

/* -------------------- Implementation: Mesher -------------------- */

inline
ResamplingMesher::ResamplingMesher (void)
{
}

inline void
ResamplingMesher::set_input_mesh (InputMeshData const& rmesh)
{
  this->rmesh = rmesh;
}

inline void
ResamplingMesher::set_config (ResamplingConf const& conf)
{
  this->config = conf;
}

inline void
ResamplingMesher::set_sampler (ResamplingSampler const& sampler)
{
  this->sampler = &sampler;
}

inline OutputMeshData
ResamplingMesher::get_output_mesh (void) const
{
  return this->emesh;
}

/* ------------------ Implementation: Resampling ------------------ */

inline
Resampling::Resampling (void)
{
}

inline void
Resampling::set_config (ResamplingConf const& config)
{
  this->config = config;
}

inline void
Resampling::set_reference_mesh (TriangleMeshPtr rmesh)
{
  this->rmesh.mesh = rmesh;
}

inline void
Resampling::set_reference_vinfo (VertexInfoListPtr rvinfo)
{
  this->rmesh.vinfo = rvinfo;
}

inline void
Resampling::set_reference_features (FeatureEdgesPtr rfeatures)
{
  this->rmesh.features = rfeatures;
}

inline void
Resampling::set_reference_density (DensityFieldPtr rdensity)
{
  this->rmesh.density = rdensity;
}

inline TriangleMeshPtr
Resampling::get_resampled_mesh (void) const
{
  return this->emesh.mesh;
}

inline FeatureEdgesPtr
Resampling::get_resampled_features (void) const
{
  return this->emesh.features;
}

inline VertexRefListPtr
Resampling::get_vertex_references (void) const
{
  return this->emesh.vrefs;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_RESAMPLING_HEADER */
