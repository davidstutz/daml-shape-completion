#include <iostream>

#include "exception.h"
#include "triangle2.h"
#include "straightline2.h"
#include "microdelaunay.h"
#include "densitytriangle2.h"

REMESHER_NAMESPACE_BEGIN

bool
MicroDelaunay::is_delaunay (void) const
{
  Vec2f center(0.0f, 0.0f);

  /* Take two subsequent vertices, build triangle with center and check
   * delaunay with next and previous vertex in the loop. */
  std::size_t as = this->adj2d.size();
  for (std::size_t i = 0; i < as; ++i)
  {
    std::size_t ip1 = (i + 1) % as;
    std::size_t ip2 = (i + 2) % as;
    std::size_t ip3 = (i + 3) % as;

    Vec2f const& a(this->adj2d[i]);
    Vec2f const& b(this->adj2d[ip1]);
    Vec2f const& c(this->adj2d[ip2]);
    Vec2f const& d(this->adj2d[ip3]);

    float circum_qradius;
    Vec2f circum_center;
    Triangle2f t(center, b, c);
    t.get_circumcircle(circum_qradius, circum_center);

    if ((a - circum_center).qlength() < circum_qradius
        || (d - circum_center).qlength() < circum_qradius)
      return false;
  }

  return true;
}

/* ---------------------------------------------------------------- */

/*
 * The area of the polygon (the voronoi cell) is calculated as
 *
 *   A = 1/2 * Sum(xi * yi+1 - xi+1 * y1)
 *
 * The centroid of the polygon is calcuated as:
 *
 *   x = 1 / 6A * Sum( (xi + xi+1) * (xi * yi+1 - xi+1 * y1) )
 *   y = 1 / 6A * Sum( (yi + yi+1) * (xi * yi+1 - xi+1 * y1) )
 */

Vec2f
MicroDelaunay::get_voronoi_centroid (void) const
{
  Polygon2 poly;
  PolygonDensity poly_dens;
  this->fill_polygon_values(poly, poly_dens);

#if 0
  /* First create the Voronoi cell polygon. */
  Polygon2 poly;
  std::size_t as = this->adj2d.size();
  for (std::size_t i = 0; i < as; ++i)
  {
    std::size_t ip1 = (i + 1) % as;
    Triangle2f t(center, this->adj2d[i], this->adj2d[ip1]);
    poly.push_back(t.get_circumcenter());
  }

  /* If there are feature edges set, the delaunay triangulation is
   * constrained. To get good results near features, the Voronoi
   * polygon is clipped with the feature edges. */
  if (!this->features.empty())
  {
    for (std::size_t i = 0; i < this->features.size(); ++i)
      if (this->features[i])
      {
        //if (verbose)
        //  std::cout << "Feature found, clipping polygon" << std::endl;
        std::size_t ip1 = (i + 1) % as;
        poly.clip_with_line(this->adj2d[i], this->adj2d[ip1]);
      }
  }

  /* If there are density values for the micro patch vertices, compute
   * density values for the Voronoi polygon vertices by locating the fan
   * area the polygon vertex belongs to, and calculating the density value
   * with linear interpolation. */
  PolygonDensity poly_dens;
  if (!this->adj_density.empty())
  {
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
      /* Find a triangle the polygon vertex belongs to. To do that
       * we use two lines along the face edges and check if the
       * vertex is left-of the right line and right-of the left line. */
      std::size_t tri = MAX_SIZE_T;
      for (std::size_t j = 0; tri == MAX_SIZE_T && j < as; ++j)
      {
        std::size_t jp1 = (j + 1) % as;
        StraightLine2f line1(center, this->adj2d[j]);
        StraightLine2f line2(center, this->adj2d[jp1]);

        if (line1.edge_equation(poly[i]) <= 0.0f
            && line2.edge_equation(poly[i]) >= 0.0f)
          tri = j;
      }

      if (tri == MAX_SIZE_T)
        throw Exception("Cannot find corresponding triangle");

      /* Get barycentric coordinate of the voronoi vertex. */
      std::size_t trip1 = (tri + 1) % as;
      Triangle2f t2(center, this->adj2d[tri], this->adj2d[trip1]);
      Vec3f bary = t2.get_bary_coords(poly[i]);

      /* Now calculate the density value for the polygon vertex. If
       * the polygon vertex is inside a triangle, use barycentric
       * coordinates for linear density interpolation. Otherwise
       * calculate the distance to the edge lines and use that distance
       * as linear interpolation weight for the fan density vertices. */
      float density;
      if (bary[0] < 0.0f || bary[1] < 0.0f || bary[2] < 0.0f)
      {
        /* The voronoi vertex is outside the patch. Approximate the density
         * for the voronoi vertex using the distance to the lines. */
        StraightLine2f line1(center, this->adj2d[tri]);
        StraightLine2f line2(center, this->adj2d[trip1]);
        float dist1 = std::sqrt(line1.point_qdist(poly[i]));
        float dist2 = std::sqrt(line2.point_qdist(poly[i]));
        float w1 = dist2 / (dist1 + dist2);
        float w2 = dist1 / (dist1 + dist2);
        density = w1 * this->adj_density[tri] + w2 * this->adj_density[trip1];
      }
      else
      {
        /* Voronoi vertex is inside the triangle. Calculation of the
         * density for the voronoi vertex is easy in this case. */
        density = this->center_density * bary[0]
            + this->adj_density[tri] * bary[1]
            + this->adj_density[trip1] * bary[2];
      }

      /* The density is powered now, and this is a emirical value
       * to get a similar distribution to the simplification. */
      poly_dens.push_back(std::pow(density, 3.0f));
      //poly_dens.push_back(density + 0.5f);
    }
  }
#endif

  /* Actually calculate the centroid now using the polygon and, if
   * specified, the density values at the polygon vertices. */
  Vec2f centroid(0.0f, 0.0f);
  if (poly_dens.empty())
  {
    /* Create the centroid of the polygon without any density
     * information, that is a uniform density distributed over
     * the polygon. The calculation is simple in this case. */
    float area = 0.0f;
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
      std::size_t ip1 = (i + 1) % poly.size();
      float part_area = poly[i][0] * poly[ip1][1] - poly[ip1][0] * poly[i][1];
      centroid[0] += (poly[i][0] + poly[ip1][0]) * part_area;
      centroid[1] += (poly[i][1] + poly[ip1][1]) * part_area;
      area += part_area;
    }

    area /= 2.0f;
    centroid /= area * 6.0f;
  }
  else
  {
    /* We have a polygon and density values for each polygon vertex.
     * We aim at calculating the centroid of the polygon. */

#if 1
    /* We calculate the mass and the centroid for each triangle
     * in the micro patch. We sum the centroids and weight them
     * with their corresponding mass. */
    float total_mass = 0.0f;
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
      std::size_t ip1 = (i + 1) % poly.size();
      Vec3f dens(this->center_density, poly_dens[i], poly_dens[ip1]);
      DensityTriangle2f tri(Vec2f(0.0f, 0.0f), poly[i], poly[ip1], dens);
      float mass = tri.get_mass();
      Vec2f c = tri.get_centroid();

      centroid += c * mass;
      total_mass += mass;
    }
    centroid /= total_mass;

#else

    /* This implementation is an approximation only FOR NOW. */
    float total_dens = 0.0f;
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
      centroid += poly[i] * poly_dens[i];
      total_dens += poly_dens[i];
    }
    centroid /= total_dens;
#endif
  }

  return centroid;
}

/* ---------------------------------------------------------------- */

Vec2f
MicroDelaunay::get_area_optimal_center (void) const
{
  MicroPatch::AdjacentVertices2D const& v = this->adj2d;
  Vec2f center(0.0f, 0.0f);

  /* Calculate total area of the polygon. */
  float parea = 0.0f;
  for (std::size_t i = 0; i < v.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % v.size();
    Triangle2f tri(center, v[i], v[ip1]);
    parea += tri.get_area();
  }

  /* Calculate the desired area for each triangle. */
  std::vector<float> tarea;
  if (this->adj_density.empty())
  {
    tarea.resize(v.size(), parea / (float)v.size());
  }
  else
  {
    tarea.resize(v.size());
    float total_density = 0.0f;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      std::size_t ip1 = (i + 1) % v.size();
      tarea[i] = 1 / (this->adj_density[i] + this->adj_density[ip1]);
      total_density += tarea[i];
    }
    for (std::size_t i = 0; i < v.size(); ++i)
      tarea[i] = parea * tarea[i] / total_density;
  }

  #if 0
  std::cout << "Weights: ";
  for (std::size_t i = 0; i < v.size(); ++i)
    std::cout << tarea[i] << ", ";
  std::cout << std::endl;
  #endif

  Matrix2f m;
  m[0] = 0.0f;
  m[1] = 0.0f;
  m[2] = 0.0f;
  m[3] = 0.0f;
  Vec2f b(0.0f, 0.0f);

  for (std::size_t i = 0; i < v.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % v.size();
    float xip1xi = (v[ip1][0] - v[i][0]);
    float yiyip1 = (v[i][1] - v[ip1][1]);

    m[0] += yiyip1 * yiyip1;
    m[1] += xip1xi * yiyip1;
    m[2] += xip1xi * yiyip1;
    m[3] += xip1xi * xip1xi;

    float xiyip1xip1yi = v[i][0] * v[ip1][1] - v[ip1][0] * v[i][1] - tarea[i];
    b[0] -= xiyip1xip1yi * yiyip1;
    b[1] -= xiyip1xip1yi * xip1xi;
  }

  float det = m[0] * m[3] - m[1] * m[2];

  //std::cout << "Linear system Ax=b: " << m << " * x = " << b << std::endl;
  //std::cout << "Determinant is " << det << std::endl;

  //if (FLOAT_EQ(det, 0.0f))
  if (det == 0.0f)
    throw Exception("Area determinant is zero!");

  Matrix2f invm;
  invm[0] = m[3] / det;
  invm[1] = -m[1] / det;
  invm[2] = -m[2] / det;
  invm[3] = m[0] / det;

  //std::cout << "Solution is " << (invm * b) << std::endl;

  return invm * b;
}

/* ---------------------------------------------------------------- */

void
MicroDelaunay::get_voronoi_area_and_mass (float& area, float& mass)
{
  Polygon2 poly;
  PolygonDensity poly_dens;
  this->fill_polygon_values(poly, poly_dens);

  if (poly_dens.empty())
  {
    poly_dens.resize(poly.size(), 1.0f);
    this->center_density = 1.0f;
  }

  area = 0.0f;
  mass = 0.0f;
  for (std::size_t i = 0; i < poly.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % poly.size();
    Vec3f dens(this->center_density, poly_dens[i], poly_dens[ip1]);
    DensityTriangle2f tri(Vec2f(0.0f, 0.0f), poly[i], poly[ip1], dens);
    mass += tri.get_mass();
    area +=  tri.get_area();
  }
}

/* ---------------------------------------------------------------- */

void
MicroDelaunay::debug_bary_coords (Vec2f const& pnt) const
{
  for (std::size_t i = 0; i < this->adj2d.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % this->adj2d.size();
    Triangle2f tri(Vec2f(0.0f, 0.0f), this->adj2d[i], this->adj2d[ip1]);
    Vec3f bc = tri.get_bary_coords(pnt);
    std::cout << "Barycentric coordinates: " << bc << std::endl;
    if (bc[0] > 0.0f && bc[1] > 0.0f && bc[2] > 0.0f)
      std::cout << "  Real coords: " << (Vec2f(0.0f, 0.0f) * bc[0]
          + this->adj2d[i] * bc[1] + this->adj2d[ip1] * bc[2]) << std::endl;
  }
}

/* ---------------------------------------------------------------- */

std::size_t
MicroDelaunay::get_face_and_bary (Vec2f const& pnt, Vec2f& bary) const
{
  for (std::size_t i = 0; i < this->adj2d.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % this->adj2d.size();
    Triangle2f tri(Vec2f(0.0f, 0.0f), this->adj2d[i], this->adj2d[ip1]);
    Vec3f bc = tri.get_bary_coords(pnt);
    if (bc[0] >= 0.0f && bc[1] >= 0.0f && bc[2] >= 0.0f)
    {
      bary = Vec2f(bc[0], bc[1]);
      return i;
    }
  }

  throw Exception("Point not inside a triangle!");
}

/* ---------------------------------------------------------------- */

void
MicroDelaunay::fill_polygon_values (Polygon2& poly, PolygonDensity& poly_dens) const
{
  poly.clear();
  poly_dens.clear();

  Vec2f center(0.0f, 0.0f); // Short hand

  /* First create the Voronoi cell polygon. */
  std::size_t as = this->adj2d.size();
  for (std::size_t i = 0; i < as; ++i)
  {
    std::size_t ip1 = (i + 1) % as;
    Triangle2f t(center, this->adj2d[i], this->adj2d[ip1]);
    poly.push_back(t.get_circumcenter());
  }

  /* If there are feature edges set, the delaunay triangulation is
   * constrained. To get good results near features, the Voronoi
   * polygon is clipped with the feature edges. */
  if (!this->features.empty())
  {
    for (std::size_t i = 0; i < this->features.size(); ++i)
      if (this->features[i])
      {
        //if (verbose)
        //  std::cout << "Feature found, clipping polygon" << std::endl;
        std::size_t ip1 = (i + 1) % as;
        poly.clip_with_line(this->adj2d[i], this->adj2d[ip1]);
      }
  }

  /* If there are density values for the micro patch vertices, compute
   * density values for the Voronoi polygon vertices by locating the fan
   * area the polygon vertex belongs to, and calculating the density value
   * with linear interpolation. */
  if (!this->adj_density.empty())
  {
    for (std::size_t i = 0; i < poly.size(); ++i)
    {
      /* Find a triangle the polygon vertex belongs to. To do that
       * we use two lines along the face edges and check if the
       * vertex is left-of the right line and right-of the left line. */
      std::size_t tri = MAX_SIZE_T;
      for (std::size_t j = 0; tri == MAX_SIZE_T && j < as; ++j)
      {
        std::size_t jp1 = (j + 1) % as;
        StraightLine2f line1(center, this->adj2d[j]);
        StraightLine2f line2(center, this->adj2d[jp1]);

        if (line1.edge_equation(poly[i]) <= 0.0f
            && line2.edge_equation(poly[i]) >= 0.0f)
          tri = j;
      }

      if (tri == MAX_SIZE_T)
        throw Exception("Cannot find corresponding triangle");

      /* Get barycentric coordinate of the voronoi vertex. */
      std::size_t trip1 = (tri + 1) % as;
      Triangle2f t2(center, this->adj2d[tri], this->adj2d[trip1]);
      Vec3f bary = t2.get_bary_coords(poly[i]);

      /* Now calculate the density value for the polygon vertex. If
       * the polygon vertex is inside a triangle, use barycentric
       * coordinates for linear density interpolation. Otherwise
       * calculate the distance to the edge lines and use that distance
       * as linear interpolation weight for the fan density vertices. */
      float density;
      if (bary[0] < 0.0f || bary[1] < 0.0f || bary[2] < 0.0f)
      {
        /* The voronoi vertex is outside the patch. Approximate the density
         * for the voronoi vertex using the distance to the lines. */
        StraightLine2f line1(center, this->adj2d[tri]);
        StraightLine2f line2(center, this->adj2d[trip1]);
        float dist1 = std::sqrt(line1.point_qdist(poly[i]));
        float dist2 = std::sqrt(line2.point_qdist(poly[i]));
        float w1 = dist2 / (dist1 + dist2);
        float w2 = dist1 / (dist1 + dist2);
        density = w1 * this->adj_density[tri] + w2 * this->adj_density[trip1];
      }
      else
      {
        /* Voronoi vertex is inside the triangle. Calculation of the
         * density for the voronoi vertex is easy in this case. */
        density = this->center_density * bary[0]
            + this->adj_density[tri] * bary[1]
            + this->adj_density[trip1] * bary[2];
      }

      /* The density is powered now, and this is a emirical value
       * to get a similar distribution to the simplification. */
      poly_dens.push_back(std::pow(density, 3.0f));
      //poly_dens.push_back(density + 0.5f);
    }
  }
}

REMESHER_NAMESPACE_END
