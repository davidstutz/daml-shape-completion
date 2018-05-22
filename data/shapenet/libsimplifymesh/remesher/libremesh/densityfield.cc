#include <cstdlib>
#include <iostream>

#include "elapsedtimer.h"
#include "exception.h"
#include "triangle3.h"
#include "densityfield.h"

REMESHER_NAMESPACE_BEGIN

#define DENS_SCALE_MIN 0.0f
#define DENS_SCALE_MAX 2.0f

void
DensityField::calculate_density_field (void)
{
  if (this->mesh.get() == 0)
    throw Exception("DensityField: No mesh has been set");

  if (this->mesh->get_vertices().empty())
    return;

  if (this->mesh->get_vertices().size() != this->vinfo->size())
    throw Exception("DensityField: Vertex information not properly set");

  #if 0
  /* Test some density gradient over the mesh. */
  MeshVertexList const& verts = this->mesh->get_vertices();
  this->clear();
  this->resize(verts.size(), VertexDensity(1.0f));
  this->max_density = 0.0f;
  Vec3f bb_min(INFINITY, INFINITY, INFINITY);
  Vec3f bb_max(-INFINITY, -INFINITY, -INFINITY);
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    bb_min[0] = (verts[i][0] < bb_min[0] ? verts[i][0] : bb_min[0]);
    bb_min[1] = (verts[i][1] < bb_min[1] ? verts[i][1] : bb_min[1]);
    bb_min[2] = (verts[i][2] < bb_min[2] ? verts[i][2] : bb_min[2]);

    bb_max[0] = (verts[i][0] > bb_max[0] ? verts[i][0] : bb_max[0]);
    bb_max[1] = (verts[i][1] > bb_max[1] ? verts[i][1] : bb_max[1]);
    bb_max[2] = (verts[i][2] > bb_max[2] ? verts[i][2] : bb_max[2]);
  }
  Vec3f centroid((bb_min + bb_max) / 2.0f);
  float scale(std::max(bb_max[0] - bb_min[0],
      std::max(bb_max[1] - bb_min[1], bb_max[2] - bb_min[2])));

  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    VertexDensity& d(this->at(i));
    Vec3f v(verts[i]);
    v -= centroid;
    v /= scale;

    d.density = std::pow(v[1] *-1.0f + 1.0f, this->config.contrast_exp * 10.0f);
    if (d.density > this->max_density)
      this->max_density = d.density;
    std::cout << "Setting density " << d.density
        << " for vertex " << i << std::endl;
  }

  return;
  #endif

  #if 0
  // FIXME: Remove this at some time, this is a test hack
  if (this->mesh->get_vertices().size() == 8)
  {
    // Might be the box model, apply some testing
    this->clear();
    this->resize(8, VertexDensity(10.0f));
    this->max_density = 10.0f;
    
    for (std::size_t i = 0; i < 8; ++i)
    {
      if (i < 4)
        this->at(i) = VertexDensity(1.0f);
      else
        this->at(i) = VertexDensity(2.0f);
    }

    return;
  }
  #endif

  ElapsedTimer t;

  std::cout << "Calculating density values for each vertex..." << std::endl;

  /* Calculate the curvature values for each vertex. The curvature
   * values are used to calculate the density value. The density value
   * is a combination of the absolute gaussian curvature and the
   * absolute mean curvature, clamped to [0.0f, 1000.0f] and
   * scaled to range [0.0f, 2.0f]. */
  this->calculate_density();

  /* Apply a contrast function. An exponent of 1.0f has no effect
   * and the contrast function is skipped in that case. */
  if (this->config.contrast_exp != 1.0f)
    this->apply_contrast_function();

  /* Do laplacian smoothing. If either iterations is zero or the smooth
   * factor is zero, smoothing has no effect and is skipped. */
  if (this->config.smooth_iter != 0 && this->config.smooth_factor != 0.0f)
    this->smooth_density_field();

  std::cout << "Calculating density values took " << t.get_elapsed()
      << "ms" << std::endl;
}

/* ---------------------------------------------------------------- */

void
DensityField::calculate_density (void)
{
  /* Create some short hands. */
  MeshVertexList const& verts = this->mesh->get_vertices();

  /* Calculate curvatures for each vertex. */
  this->max_density = 0.0f;
  this->clear();
  this->resize(this->vinfo->size());
  for (std::size_t i = 0; i < this->vinfo->size(); ++i)
  {
    /* Request adjacent vertices. */
    VertexInfo::VertexList vlist;
    this->vinfo->get_adjacent_vertices(i, vlist);

    /* Treat vertices on features creases differently (if requested). */
    float feature_curvature = 0.0f;
    bool is_feature_vertex = false;
    if (this->features.get() && !this->features->empty())
    {
      FeatureVertexEdges fve = this->features->get_edges_for_vertex(i, vlist);

      if (fve.size() > 0)
        is_feature_vertex = true;

      if (this->config.features_by_angle
          && !this->config.no_feature_density
          && fve.size() == 2)
      {
        Vec3f e1 = (verts[i] - verts[fve[0]]).norm();
        Vec3f e2 = (verts[fve[1]] - verts[i]).norm();
        float cos_angle = e1 * e2;
        float real_angle = std::acos(cos_angle);
        float real_thres_angle = std::acos(this->config.feature_angle);
        float angle_multiplier = std::min(real_angle / real_thres_angle, 1.0f);
        feature_curvature = this->config.max_curvature * angle_multiplier;
      }
    }

    /* If feature density is unwanted, set curvature to min. */
    if (is_feature_vertex && this->config.no_feature_density)
    {
      this->at(i) = VertexDensity(this->config.min_curvature);
      continue;  
    }

    /* Leave the uninitialized density value for some vertex classes. */
    if (this->vinfo[i].vclass != VERTEX_CLASS_SIMPLE)
      continue;

    if (vlist.size() < 3)
    {
      std::cout << "Warning: Simple vertex " << i << " with "
          << this->vinfo[i].adj_faces.size()
          << " < 3 adjacent faces!" << std::endl;
      this->at(i).valid = false;
      continue;
      //throw Exception("Invalid simple vertex with <3 neighbours");
    }

    /* Start computation by iterating over the triangles. This computation
     * determines the face areas and integral gaussian and mean curvature. */
    float face_areas = 0.0f;
    float int_gaussian = 2.0f * (float)MY_PI;
    float int_abs_mean = 0.0f;

    /* Three edges with its length and two normals are maintained. */
    Vec3f e1;
    Vec3f e2 = verts[vlist[0]] - verts[i];
    Vec3f e3 = verts[vlist[1]] - verts[i];
    float len_e1;
    float len_e2 = e2.length();
    Vec3f n12;
    Vec3f n23 = e2.cross(e3).norm();

    for (std::size_t j = 0; j < vlist.size(); ++j)
    {
      /* Update edges, length and normals. */
      e1 = e2;
      e2 = e3;
      e3 = verts[vlist[(j + 2) % vlist.size()]] - verts[i];

      len_e1 = len_e2;
      len_e2 = e2.length();
      float len_e12 = (e2 - e1).length();

      n12 = n23;
      n23 = e2.cross(e3).norm();

      /* Compute area of the triangle. */
      float p = (len_e1 + len_e2 + len_e12) / 2.0f;
      float qfarea = p * (p - len_e1) * (p - len_e2) * (p - len_e12);
      qfarea = MY_MAX(0.0f, qfarea);
      face_areas += std::sqrt(qfarea);

      /* Compute integral gaussian curvaure. */
      {
        float edge_scalar = e1.scalar(e2) / (len_e1 * len_e2);
        edge_scalar = MY_MIN(1.0f, edge_scalar);
        edge_scalar = MY_MAX(-1.0f, edge_scalar);
        float angle = ::acosf(edge_scalar);
        int_gaussian -= angle;
      }

      /* Compute integral absolute mean curvature. */
      {
        float normal_scalar = n12.scalar(n23);
        normal_scalar = MY_MIN(1.0f, normal_scalar);
        normal_scalar = MY_MAX(-1.0f, normal_scalar);
        int_abs_mean += 0.25f * len_e2 * ::fabsf(::acosf(normal_scalar));
      }

#if REMESHER_NAN_CHECKS
      /* Check for NAN values to detect problems. */
      if (std::isnan(int_abs_mean) || std::isnan(int_gaussian)
          || std::isnan(face_areas))
      {
        std::cout << "DensityField: NAN in iteration " << j << std::endl
            << "  Absolute mean: " << int_abs_mean << std::endl
            << "  Gaussian: " << int_gaussian << std::endl
            << "  Face areas: " << face_areas << ", p: " << p
            << ", len1: " << len_e1 << ", len2: " << len_e2
            << ", len12: " << len_e12 << std::endl
            << "  p-l1: " << (p - len_e1) << std::endl;
        ::exit(1);
      }
#endif
    }

    if (face_areas <= 0.0f)
    {
      std::cout << "Warning: Vertex area is " << face_areas << std::endl;
      this->at(i).valid = false;
      continue;
    }

    /* Compute gaussian and absolute mean curvature. The barycentric area
     * of the triangles is used, therefore one third of the area is used. */
    float gaussian = 3.0f * int_gaussian / face_areas;
    float abs_mean = 3.0f * int_abs_mean / face_areas;

#if REMESHER_NAN_CHECKS
    if (std::isnan(abs_mean) || std::isnan(gaussian))
      std::cout << "Warning: Curvature is NAN: abs_mean: " << abs_mean
          << ", gaussian: " << gaussian << std::endl;
#endif

    float gaussian_part = this->config.alpha * std::fabs(gaussian);
    float abs_mean_part = this->config.beta * abs_mean * abs_mean;
    float curvature = gaussian_part + abs_mean_part;

    /* Clamp the curvature to the configured range. */
    curvature = MY_MAX(this->config.min_curvature, curvature);
    curvature = MY_MIN(this->config.max_curvature, curvature);

    /* Take feature curvature into consideration. */
    curvature = MY_MAX(curvature, feature_curvature);

    /* Scale to curvature value to [0.0, 1.0]. */
    float density = (curvature - this->config.min_curvature)
        / (this->config.max_curvature - this->config.min_curvature);

    /* Scale the density to [SCALE_MIN, SCALE_MAX]. */
    density = density * (DENS_SCALE_MAX - DENS_SCALE_MIN) + DENS_SCALE_MIN;

    /* Assign the density. */
    this->at(i) = VertexDensity(density);
    if (density > this->max_density)
      this->max_density = density;
  }
}

/* ---------------------------------------------------------------- */

void
DensityField::smooth_density_field (void)
{
  if (this->config.smooth_iter == 0)
    return;

  std::cout << "Smoothing density field: " << std::flush;

  for (std::size_t iter = 0; iter < this->config.smooth_iter; ++iter)
  {
    if (iter % 10 == 0 && iter != 0)
      std::cout << iter << std::flush;
    else
      std::cout << "." << std::flush;

    /* Backup density information. */
    std::vector<float> density;
    density.resize(this->size());
    for (std::size_t i = 0; i < this->size(); ++i)
      density[i] = this->at(i).density;

    /* Calculate new density. */
    for (std::size_t i = 0; i < this->size(); ++i)
    {
      /* Skip smoothing for invalid vertices. */
      if (!this->at(i).valid)
        continue;

      VertexInfo::VertexList vlist;
      this->vinfo->get_adjacent_vertices(i, vlist);

      float density_sum = 0.0f;
      std::size_t valid_neighbours = 0;
      for (std::size_t v = 0; v < vlist.size(); ++v)
      {
        /* Skip influence of invalid neighbours. */
        if (!this->at(vlist[v]).valid)
          continue;

        density_sum += density[vlist[v]] - density[i];
        valid_neighbours += 1;
      }

      /* Update with smoothed density information. */
      if (valid_neighbours > 0)
      {
        float new_dens = density[i] + this->config.smooth_factor
            * density_sum / (float)valid_neighbours;
        new_dens = MY_MAX(DENS_SCALE_MIN, new_dens);
        new_dens = MY_MIN(DENS_SCALE_MAX, new_dens);
        this->at(i).density = new_dens;
      }
    }
  }

  std::cout << " done." << std::endl;
}

/* ---------------------------------------------------------------- */

void
DensityField::apply_contrast_function (void)
{
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    if (!this->at(i).valid)
      continue;
    float density = this->at(i).density;
    this->at(i).density = std::pow(density, this->config.contrast_exp);
  }
}

/* ---------------------------------------------------------------- */

VertexDensity
DensityField::get_density (std::size_t face, Vec2f const& bary)
{
  MeshFaceList const& faces = this->mesh->get_faces();
  std::size_t v1 = faces[face * 3 + 0];
  std::size_t v2 = faces[face * 3 + 1];
  std::size_t v3 = faces[face * 3 + 2];

  if (this->at(v1).valid && this->at(v2).valid && this->at(v3).valid)
  {
    float density = this->at(v1).density * bary[0]
        + this->at(v2).density * bary[1]
        + this->at(v3).density * (1.0f - bary[0] - bary[1]);
    return VertexDensity(density);
  }

  return VertexDensity();
}

/* ---------------------------------------------------------------- */

void
DensityField::get_face_info (std::size_t face,
    float* area, float* density) const
{
  Triangle3f tri(this->mesh, face);
  float _area = tri.get_area();
  float _dens = _area;

  if (!this->empty() && density != 0)
  {
    float total_dens = 0.0f;
    for (int i = 0; i < 3; ++i)
    {
      VertexDensity const& vd(this->at(this->mesh->get_faces()[face * 3 + i]));
      total_dens += (vd.valid ? vd.density : 1.0f);
    }
    _dens = (1.0f / 3.0f) * _area * total_dens;
  }

  /* Set area and density values. */
  if (area)
    *area = _area;
  if (density)
    *density = _dens;
}

/* ---------------------------------------------------------------- */

void
DensityField::get_edge_info (std::size_t v1, std::size_t v2,
    float* length, float* density) const
{
  Vec3f const& _v1(this->mesh->get_vertices()[v1]);
  Vec3f const& _v2(this->mesh->get_vertices()[v2]);
  float _len = (_v2 - _v1).length();
  float _dens = _len;

  if (!this->empty() && density != 0)
  {
    VertexDensity const& vd1(this->at(v1));
    VertexDensity const& vd2(this->at(v2));
    float total_dens = (vd1.valid ? vd1.density : 1.0f)
        + (vd2.valid ? vd2.density : 1.0f);
    _dens = (1.0f / 2.0f) * _len * total_dens;
  }

  /* Set length and density values. */
  if (length)
    *length = _len;
  if (density)
    *density = _dens;
}

REMESHER_NAMESPACE_END
