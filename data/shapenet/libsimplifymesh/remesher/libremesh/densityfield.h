#ifndef REMESHER_DENSITY_FIELD_HEADER
#define REMESHER_DENSITY_FIELD_HEADER

#include <vector>

#include "defines.h"
#include "refptrarray.h"
#include "vec2.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

struct DensityFieldConf
{
  /* Alpha and beta describes how gauss and abs. mean curvature
   * values are weighted. Alpha and beta are greater zero and sum to
   * unity. Typical values are alpha = beta = 0.5f. */
  float alpha;
  float beta;

  /* Curvature values are clamped to range to erase very high values.
   * The minimum curvauture is typically set to 0.0f, the maximum
   * curvature is typically about 1000.0f, depending on the model. */
  float min_curvature;
  float max_curvature;

  /* The exponent for the contrast function, which is basically a gamma
   * function. An exponent of 1.0f disables the contrast function. */
  float contrast_exp;

  /* Laplacian smoothing settings. Smooth factor describes the intensity
   * of smoothing, typically around 0.5f. Smoothing is executed iteration
   * times, causing in values being distributed over more edges. */
  float smooth_factor;
  std::size_t smooth_iter;

  /* Determines if vertices lying on a crease should be treated differently.
   * Corner vertices are still regularly handled, density for inner edge
   * vertices are evaluated using the angle between the edges. The angle
   * is to be specified in cos(angle), where angle is in RAD. */
  bool features_by_angle;
  float feature_angle;

  /** Specifies if vertices on features should get min density assigned. */
  bool no_feature_density;

  DensityFieldConf (void);
};

/* ---------------------------------------------------------------- */

/*
 * A class to store the gaussian and absolute mean curvature.
 * The valid flag is set the vertex does not belong to a valid facette.
 */
struct VertexDensity
{
  // REMOVE valid? Use negative density to indicate invalid density?
  bool valid;
  //float gaussian;
  //float abs_mean;
  float density;

  VertexDensity (void);
  VertexDensity (float density);
  //VertexDensity (float gaussian, float abs_mean);
};

/* ---------------------------------------------------------------- */

/*
 * The density field calculates a curvature value for each vertex.
 * If the vertex neighbours are supplied, this data structure is
 * used. Otherwise vertex neighbours are calculated. The curvature
 * values are used to calculate a density value for the vertex,
 * see parameters "alpha", "beta", and "contrast_exp". Smoothing
 * is performed to distribute density values and to produce a
 * smooth graduation between highly curved and flat regions.
 *
 * The method to calculate the absolute mean curvature and the
 * gaussian curvature is described in:
 *
 *   Nira Dyn, Kai Hormann, Sun-Jeong Kim, and David Levin
 *   Optimizing 3D Triangulations Using Discrete Curvature Analysis
 *
 * Where the barycentric area is used for area calculation.
 *
 * TODO: Density values for boundary vertices
 */

class DensityField;
typedef RefPtrArray<DensityField, VertexDensity> DensityFieldPtr;

class DensityField : public std::vector<VertexDensity>
{
  private:
    DensityFieldConf config;
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    FeatureEdgesPtr features;
    float max_density;

  protected:
    DensityField (void);

    void calculate_density (void);
    void smooth_density_field (void);
    void apply_contrast_function (void);

  public:
    static DensityFieldPtr create (void);
    DensityFieldPtr create_copy (void) const;

    /* Set algorithm configuration. */
    void set_config (DensityFieldConf const& conf);

    /* Set algorithm data (mesh and vertex info required). */
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);
    void set_feature_edges (FeatureEdgesPtr features);

    /* Start algorithm; calculate density values. */
    void calculate_density_field (void);

    /* Returns the maxiumum density value present; value is in [0, max].
     * May be used to scale each value to e.g. a color range [0, 1]. */
    float get_max_density (void) const;

    /* Provides a density value for a point inside a face. */
    VertexDensity get_density (std::size_t face, Vec2f const& bary);

    /* Provides area and integrated density of a face of the mesh.
     * The arguments may be NULL if no result is required.
     * The density equals the area if the density field is uninitialized. */
    void get_face_info (std::size_t face, float* area, float* density) const;
    float get_face_density (std::size_t face) const;

    /* Provides length and integrated density of an edge of the mesh.
     * Arguments may be NULL. Density equals length if
     * the density field is uninitialized. */
    void get_edge_info (std::size_t v1, std::size_t v2,
        float* length, float* density) const;
    float get_edge_density (std::size_t v1, std::size_t v2) const;

    /* Returns the amount of memory used. */
    std::size_t get_memory_usage (void) const;
};

/* ================================================================ */

inline
DensityFieldConf::DensityFieldConf (void)
{
  this->alpha = 0.5f;
  this->beta = 0.5f;
  this->min_curvature = 0.0f;
  this->max_curvature = 1000.0f;
  this->contrast_exp = 1.0f;
  this->smooth_factor = 0.5f;
  this->smooth_iter = 10;
  this->features_by_angle = true;
  this->feature_angle = std::cos((float)MY_DEG2RAD(10.0f));
  this->no_feature_density = false;
}

inline
VertexDensity::VertexDensity (void)
{
  this->valid = false;
  //this->gaussian = 0.0f;
  //this->abs_mean = 0.0f;
  this->density = 0.0f;
}

inline
VertexDensity::VertexDensity (float density)
{
  this->valid = true;
  this->density = density;
}

inline
DensityField::DensityField (void)
{
}

inline DensityFieldPtr
DensityField::create (void)
{
  return DensityFieldPtr(new DensityField);
}

inline DensityFieldPtr
DensityField::create_copy (void) const
{
  return DensityFieldPtr(new DensityField(*this));
}

inline void
DensityField::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
DensityField::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
DensityField::set_feature_edges (FeatureEdgesPtr features)
{
  this->features = features;
}

inline void
DensityField::set_config (DensityFieldConf const& conf)
{
  this->config = conf;
}

inline float
DensityField::get_max_density (void) const
{
  return this->max_density;
}

inline float
DensityField::get_face_density (std::size_t face) const
{
  float density;
  this->get_face_info(face, 0, &density);
  return density;
}

inline float
DensityField::get_edge_density (std::size_t v1, std::size_t v2) const
{
  float density;
  this->get_edge_info(v1, v2, 0, &density);
  return density;
}

inline std::size_t
DensityField::get_memory_usage (void) const
{
  return this->capacity() * sizeof(VertexDensity);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_DENSITY_FIELD_HEADER */
