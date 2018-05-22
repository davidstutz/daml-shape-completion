#ifndef SIMPLIFICATION_HEADER
#define SIMPLIFICATION_HEADER

#include "defines.h"
#include "vertexinfo.h"
#include "plane.h"
#include "densityfield.h"
#include "trianglemesh.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

/*
 * The simiplification in this class implements the method described in:
 *
 *   William J. Schroeder, Jonathan A. Zarge, William E. Lorensen
 *   Decimation of Triangle Meshes
 *
 * The simplification iterates in several passes over the remaining
 * vertices and removes all vertices that pass all of the performed checks.
 * With each iteration, the checks are relaxed until the desired amount
 * of vertices is reached. Remaining vertices are not moved and can
 * be associated with vertices in the original mesh.
 *
 * Besides the fidelity checks described in this paper, new checks
 * has been implemented regarding the area associated with each vertex.
 * With this addition vertex areas are bounded, wich produces a nice
 * initial vertex partition for isotropic remeshing algorithms.
 */

/* ---------------------------------------------------------------- */

struct SimplificationConf
{
  /* This is the stop criterion, the desired amount of vertices. */
  std::size_t vertex_limit;

  /* The first pass starts with the initial threshold, the threshold
   * of each subsequent pass is multiplied by threshold factor up to
   * maximal threshold. Several passes ensure uniform simplification. */
  float initial_threshold;
  float maximum_threshold;
  float threshold_factor;

  /* Defines if checks are done to ensure fidelity of the output mesh.
   * This alone typically produces non-uniform vertex distribution. */
  bool perform_fidelity_checks;
  float fidelity_check_penalty;

  /* Defines if checks are done regarding the area of triangles, thus
   * vertices with small adjacent triangles are removed first.
   * This produces uniform vertex distribution. */
  bool perform_area_checks;
  float area_check_penalty;

  /* Defines if per-vertex checks for sharp surrounding edges are done.
   * This efficently allows edge-to-vertex distance checks instead of
   * plane-to-vertex distance checks for inner edge vertices. Corner
   * vertices are still evaluated with the plane-to-vertex distance.
   * This is a useful setting even for meshes without global features.
   * The local feature angle is to be specified as cos(angle). */
  bool check_local_features;
  float local_feature_angle;

  /* Prevents vertices on global feature creases to be simplified.
   * This is typically not very useful since features should also
   * be simplified to get a proper sampling. */
  bool keep_global_features;

  SimplificationConf (void);
};

/* ---------------------------------------------------------------- */

class Simplification
{
  private:
    typedef std::vector<bool> VertexDeletedList;
    typedef std::pair<std::size_t, std::size_t> SplitLine;

  public:
    typedef std::vector<std::size_t> VertexIndexRelocList;

  private:
    SimplificationConf config;

    TriangleMeshPtr mesh;
    DensityFieldPtr vdens;
    FeatureEdgesPtr features;

    VertexInfoListPtr vinfo;
    VertexDeletedList dinfo;
    VertexIndexRelocList* ireloc;

    /* Vertex amount need to be tracked since it's not
     * available from the triangle mesh during simplification. */
    std::size_t vertex_amount;

  protected:
    /* Info preparation consists of creating a data structure provided
     * by the VertexInfoList class, classifying all vertices and
     * creating the additional vertex info for simplification. */
    void prepare_vertex_info (void);

    /* Runs the simplification process. The process has several passes. */
    void run_simplification (void);

    /* A single pass of the simplification process using threshold
     * for distance to edge, distance to plane or area checks.
     * The method returns true if there are no more vertex candidates. */
    bool run_pass (float threshold);

    /* Simplification strategies for each vertex type. */
    bool handle_removal_candidate (std::size_t index, float thres);

    /* Checks for global features for the vertex. The function returns
     * false if global features forbid the deletion of the vertex.
     * If a possible simplification has been found (for inner edge
     * vertices) the splitline is specified. */
    bool check_global_features (std::size_t index, SplitLine& sline,
        VertexInfo::VertexList const& vlist);

    /* Performs checks if simplification harms fidelity. The check
     * return true if the vertex is good to remove. The method also
     * returns the average for adjacent vertices. */
    bool check_candidate_fidelity (std::size_t index, float thres,
        VertexInfo::VertexList const& vlist, Plane3f const& plane,
        SplitLine& sline);

    /* Performs checks to ensure that only small trinalges are removed.
     * Returns true if the vertex is surrounded by resonable small
     * triangles. */
    bool check_candidate_area (std::size_t index, float thres);

    /* Replaces a set of faces (aflist) with new faces given
     * by a list of vertices in groups of three vertices each face. */
    void replace_faces (VertexInfo::FaceList aflist, MeshFaceList const& flist);

    /* Efficiently cleans removed vertices and faces from the mesh. */
    void clean_mesh (void);

    /* Debugging. */
    void memory_debug (void);

  public:
    Simplification (void);

    /* Sets the simplification configuration. If no configuration
     * is set, a default configuration is used. */
    void set_config (SimplificationConf const& conf);

    /* Sets the triangle mesh to be simplified. Simplification
     * is done in-place and will change the mesh if copy mesh
     * is not set. */
    void set_mesh (TriangleMeshPtr mesh, bool copy_mesh = false);

    /* Sets the density field for density-based simplification.
     * If the density field is omitted, uniform density is assumed. */
    void set_density_field (DensityFieldPtr density);

    /* Sets a set of feature edges. During simplification the set of
     * feature edges is modified in-place and the passed data will be
     * changed if it is not copied. */
    void set_features (FeatureEdgesPtr features, bool copy_features = false);

    /* Allows to request a vertex index relocation list. The
     * resulting list contains the indices of the vertices
     * prior simiplification. This is useful to access information
     * on the original mesh. */
    void fill_vertex_reloc_list (VertexIndexRelocList& irlist);

    /* Create the data structure, runs the simplification and cleans the
     * resulting mesh from deleted vertices. */
    void start_simplification (void);

    /* The resulting mesh can be obtained with this method. */
    TriangleMeshPtr get_mesh (void);

    /* The resulting set of feature edges can also be obtained. */
    FeatureEdgesPtr get_features (void);
};

/* ---------------------------------------------------------------- */

inline
SimplificationConf::SimplificationConf (void)
{
  this->vertex_limit = 1000;

  this->initial_threshold = 0.00001f;
  this->threshold_factor = 1.4f;
  this->maximum_threshold = 1.0f;

  this->perform_fidelity_checks = true;
  this->fidelity_check_penalty = 1.0f;
  this->perform_area_checks = true;
  this->area_check_penalty = 1.0f;

  this->check_local_features = true;
  this->local_feature_angle = 0.75f;

  this->keep_global_features = false;
}

inline void
Simplification::set_config (SimplificationConf const& conf)
{
  this->config = conf;
}

inline void
Simplification::set_density_field (DensityFieldPtr density)
{
  this->vdens = density;
}

inline void
Simplification::fill_vertex_reloc_list (VertexIndexRelocList& irlist)
{
  this->ireloc = &irlist;
}

inline TriangleMeshPtr
Simplification::get_mesh (void)
{
  return this->mesh;
}

inline FeatureEdgesPtr
Simplification::get_features (void)
{
  return this->features;
}

REMESHER_NAMESPACE_END

#endif /* SIMPLIFICATION_HEADER */
