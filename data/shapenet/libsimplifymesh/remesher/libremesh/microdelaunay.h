#ifndef REMESHER_MICRO_DELAUNAY_HEADER
#define REMESHER_MICRO_DELAUNAY_HEADER

#include "defines.h"
#include "polygon2.h"
#include "micropatch.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class takes a center vertex and surrounding vertices in 2D.
 * It calculates the voronoi centroid with or without guidance by
 * a density function and feature edges. The density values need
 * to be specified for each vertex in the micro patch.
 */
class MicroDelaunay : public MicroPatch
{
  private:
    typedef std::vector<float> AdjacentDensity;
    typedef std::vector<bool> PolygonEdges;
    typedef std::vector<float> PolygonDensity;

  private:
    float center_density;
    AdjacentDensity adj_density;
    PolygonEdges features;

  private:
    void fill_polygon_values (Polygon2& poly, PolygonDensity& poly_dens) const;

  public:
    void set_center_density (float density);
    void append_adjacent_density (float density);
    void set_feature_edge (std::size_t index);
    void clear_adjacent_density (void);

    /* Makes a sanity check if the pach is delaunay. Since due to
     * distortion the patch might not be delaunay, don't use it! */
    bool is_delaunay (void) const;

    /* Returns the Voronoi centroid, that is for Lloyd relaxation. */
    Vec2f get_voronoi_centroid (void) const;

    /* Returns area based optimal center vertex for Area Based remeshing. */
    Vec2f get_area_optimal_center (void) const;

    /* Returns the area and the mass of the Voronoi cell. */
    void get_voronoi_area_and_mass (float& area, float& mass);

    /* Returns the face ID and fills bary with the barycentric coordinate. */
    std::size_t get_face_and_bary (Vec2f const& pnt, Vec2f& bary) const;

    void debug_bary_coords (Vec2f const& pnt) const;
};

/* ---------------------------------------------------------------- */

inline void
MicroDelaunay::set_center_density (float density)
{
  this->center_density = density;
}

inline void
MicroDelaunay::append_adjacent_density (float density)
{
  this->adj_density.push_back(density);
}

inline void
MicroDelaunay::set_feature_edge (std::size_t index)
{
  if (features.empty())
    features.resize(this->adj3d.size(), false);
  features[index] = true;
}

inline void
MicroDelaunay::clear_adjacent_density (void)
{
  this->adj_density.clear();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MICRO_DELAUNAY_HEADER */
