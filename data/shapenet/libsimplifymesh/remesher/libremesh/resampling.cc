#include <cmath>
#include <iostream>

#include "helpers.h"
#include "exception.h"
#include "triangle3.h"
#include "meshskeleton.h"
#include "incdelaunay.h"
#include "delaunayflips.h"
#include "gslpolysolve.h"
#include "resampling.h"

REMESHER_NAMESPACE_BEGIN

/*
 * TODO
 * - Config option to specify feature and surface samples separately
 * - Adaptive triangle sampling
 */

void
ResamplingSampler::exec_resampling (std::size_t sample_amount)
{
  // std::cout << "Resampling: Calibrating the mesh sampler..." << std::endl;

  /* Build the mesh skeleton for feature processing. */
  this->skel = MeshSkeleton::create();
  this->skel->set_mesh(this->rmesh.mesh);
  this->skel->set_feature_edges(this->rmesh.features);
  this->skel->set_vertex_info(this->rmesh.vinfo);
  //this->skel->set_max_angle((float)MY_DEG2RAD(45.0f));
  this->skel->extract();
  this->skel->debug_print();

  /* Calculate the total area (uniform) or mass (adapted) for faces. */
  float surface_area, surface_density;
  this->calc_surface_density(surface_area, surface_density);
  float surface_scale = surface_density / surface_area;

  /* Calculate total area/density for features. */
  float feature_length, feature_density;
  this->calc_feature_density(feature_length, feature_density);
  float feature_scale = feature_density / feature_length;
  //feature_scale /= 2.0f;

  /* Calculate the number of corner vertices. */
  std::size_t corners = this->skel->get_corner_amount();

  /* Mesh density and feature stats. */
  // std::cout << "The mesh has total density "
      // << surface_density << " Ds + "
      // << feature_density << " Df, "
      // << corners << " corners, "
      // << this->skel->size() << " backbones." << std::endl;

  /* Calc vertex amount for the surface and the features. */
  /* This is done by solving the abc-formula for ax^2 + bx + c = 0. */
  std::size_t feature_samples = 0;
  std::size_t surface_samples = 0;

  if (feature_density > 0.0f)
  {
    float sqrt3 = std::sqrt(3.0f);

    float dfds = feature_scale * feature_scale / surface_scale;
    float a = 2.0f / sqrt3 * surface_density * dfds;
    float b = feature_density;
    float c = (float)corners - (float)sample_amount;

    float rf = (-b + std::sqrt(b * b - 4.0f * a * c)) / (2.0f * a);
    float rfdf = rf * feature_density;

    feature_samples = (std::size_t)MY_FLT_ROUND(rfdf);
    surface_samples = sample_amount - corners - feature_samples;
  }
  else
  {
    feature_samples = 0;
    surface_samples = sample_amount;
  }

  // std::cout << "Surface samples: " << surface_samples
  //     << ", feature samples: " << feature_samples
  //     << ", corners: " << corners << ", total samples: "
  //     << surface_samples + feature_samples + corners
  //     << std::endl;

  /* Sample corners. */
  this->sample_corners();
  /* Resample all feature creases. */
  std::size_t left = this->sample_features(feature_samples, feature_density);
  /* Release feature skeleton data. */
  this->skel.reset();

  /* Prepare algorithm output. Per edge info is created on demand. */
  this->per_face.resize(this->rmesh.mesh->get_faces().size() / 3);
  /* Resample the mesh smooth surface. */
  this->sample_surface(surface_samples + left, surface_density);
}

/* ---------------------------------------------------------------- */

void
ResamplingSampler::calc_surface_density (float& area, float& density) const
{
  MeshFaceList const& faces = this->rmesh.mesh->get_faces();
  std::size_t face_amount = faces.size() / 3;

  area = 0.0f;
  density = 0.0f;

  /* Iterate over all faces and sum area/mass. */
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    float _area, _dens;
    this->rmesh.density->get_face_info(i, &_area, &_dens);
    area += _area;
    density += _dens;
    //density += this->get_face_density(i);
  }
}

/* ---------------------------------------------------------------- */

void
ResamplingSampler::calc_feature_density (float& length, float& density) const
{
  length = 0.0f;
  density = 0.0f;

  if (this->rmesh.features->empty())
    return;

  for (std::size_t i = 0; i < this->rmesh.features->size(); ++i)
  {
    FeatureVertexEdges const& edges = this->rmesh.features[i];
    for (std::size_t j = 0; j < edges.size(); ++j)
    {
      float _len, _dens;
      this->rmesh.density->get_edge_info(i, edges[j], &_len, &_dens);
      length += _len;
      density += _dens;
      //density += this->get_edge_density(i, edges[j]);
    }
  }
}

/* ---------------------------------------------------------------- */

std::size_t
ResamplingSampler::sample_corners (void)
{
  /* Collect corner vertices. */
  std::set<std::size_t> corners;
  for (std::size_t i = 0; i < this->skel->size(); ++i)
  {
    if (this->skel->at(i).closed)
      continue;
    corners.insert(this->skel->at(i).verts.front());
    corners.insert(this->skel->at(i).verts.back());
  }

  /* Insert corner vertices to the output data structure. */
  std::set<std::size_t>::iterator iter;
  for (iter = corners.begin(); iter != corners.end(); ++iter)
    this->per_edge.push_back(SamplesPerEdge(*iter, *iter, 0));

  return corners.size();
}

/* ---------------------------------------------------------------- */

std::size_t
ResamplingSampler::sample_features (std::size_t samples, float total_density)
{
  /* Track remaining samples. */
  std::size_t remaining_samples = samples;

  /* Calculate density spacing as follows. With O the amount of
   * open backbones, S the amount of samples, D the total density
   * and C the amount of corners, the density spacing DS is calcualted
   * as follows: DS = D / (O + (S - C))
   * Corners are first sampled (and already removed from the budget) (S - C).
   */
  float dspacing = total_density / ((float)this->skel->get_open_amount()
      + (float)remaining_samples);

  // std::cout << "Sampling features with " << samples
  //     << " samples onto " << total_density << " total density "
  //     << " (spacing " << dspacing << ")" << std::endl;

  /* Sample each backbone of the mesh skeleton. */
  float error_teleport = 0.0f;
  for (std::size_t i = 0; i < skel->size(); ++i)
  {
    MeshBackbone const& mb = skel[i];

    /* First, walk over the backbone and measure density. */
    float backbone_density = this->get_backbone_density(mb);

    /* Calculate the amount of samples to distribute and the error teleport. */
    float spend_density = backbone_density + error_teleport;
    std::size_t backbone_samples = 0;
    if (mb.closed)
    {
      float flt_samples = spend_density / dspacing;
      /* Open backbones require at least three samples. Backbones
       * with less samples are considered unimportant and are ignored. */
      if (flt_samples < 2.5f)
      {
        error_teleport = spend_density;
        continue;
      }

      backbone_samples = (std::size_t)MY_FLT_ROUND(flt_samples);
      error_teleport = spend_density - (float)backbone_samples * dspacing;
    }
    else
    {
      float flt_samples = spend_density / dspacing - 1.0f;
      if (flt_samples < 0.5f)
      {
        error_teleport = spend_density;
        continue;
      }

      backbone_samples = (std::size_t)MY_FLT_ROUND(flt_samples);
      error_teleport = spend_density - (float)(backbone_samples + 1) * dspacing;
    }

    /* Finally sample the backbone. */
    // Warning: this may break ">=3 verts for closed backbones" constrained
    std::size_t distribute = std::min(backbone_samples, remaining_samples);
    if (distribute != backbone_samples)
    {
      // std::cout << "WARNING: Desired samples: " << backbone_samples
      //     << ", remaining samples: " << remaining_samples << std::endl;
    }

    this->sample_backbone(mb, backbone_density, distribute);
    remaining_samples -= distribute;
  }

  // std::cout << "Sampled feature creases with " << samples - remaining_samples
  //     << " samples, " << remaining_samples << " of " << samples
  //     << " remaining." << std::endl;

  return remaining_samples;
}

/* ---------------------------------------------------------------- */

std::size_t
ResamplingSampler::sample_surface (std::size_t samples, float total_density)
{
  MeshFaceList const& faces = this->rmesh.mesh->get_faces();
  std::size_t face_amount = faces.size() / 3;

  // std::cout << "Sampling the surface with " << samples
  //     << " samples onto " << total_density << " total density" << std::endl;

  float dens_per_sample = total_density / (float)samples;

  float teleport = 0.0f;
  std::size_t sample_amount = 0;
  float density_spent = 0.0f;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    /* Triangle density. */
    float tdens = this->rmesh.density->get_face_density(i);

    /* Sampling density is the triangle density plus teleport. */
    float sdens = tdens + teleport;
    float fsamples = sdens / dens_per_sample;
    std::size_t tsamples = (std::size_t)MY_FLT_ROUND(fsamples);

    /* Store the amount of samples in the output data. */
    this->per_face[i] = tsamples;
    sample_amount += tsamples;

    /* Calculate new teleport. */
    teleport = sdens - (float)tsamples * dens_per_sample;
    density_spent += (float)tsamples * dens_per_sample;
  }

  // std::cout << "Sampled surface with a total of " << sample_amount
  //     << " vertices on " << density_spent << " density" << std::endl;

  return 0;
}

/* ---------------------------------------------------------------- */

void
ResamplingSampler::sample_backbone (MeshBackbone const& mb,
    float backbone_dens, std::size_t samples)
{
  // std::cout << "Sampling " << (mb.closed ? "closed" : "open")
  //     << " backbone " << mb.verts.front() << "->" << mb.verts.back()
  //     << " with " << samples << " samples..." << std::endl;

  /* Hack to increase samples for closed backbones. */
  //if (mb.closed) samples *= 2;

  /* Calculate optimal density spacing for the backbone. */
  float density_spacing;
  if (mb.closed)
    density_spacing = backbone_dens / (float)samples;
  else
    density_spacing = backbone_dens / (float)(samples + 1);

  /* Set initial excess teleport. */
  float excess_teleport = 0.0f;
  if (!mb.closed)
    excess_teleport = density_spacing;

  for (BackboneIter i = mb.verts.begin(); i != mb.verts.end(); ++i)
  {
    BackboneIter ip1 = ++BackboneIter(i);
    if (ip1 == mb.verts.end())
    {
      if (!mb.closed)
        break;
      ip1 = mb.verts.begin();
    }

    std::size_t sampled = this->sample_edge(*i, *ip1, density_spacing,
        excess_teleport, samples);

    samples -= sampled;
    if (samples == 0)
      return;    
  }
}

/* ---------------------------------------------------------------- */

std::size_t
ResamplingSampler::sample_edge (std::size_t v1, std::size_t v2,
    float density_spacing, float& excess_teleport, std::size_t max_vertices)
{
  MeshVertexList const& rverts = this->rmesh.mesh->get_vertices();

  /* Length, start and end density for the edge. */
  float len = (rverts[v1] - rverts[v2]).length();
  float d1 = 1.0f;
  float d2 = 1.0f;
  if (!this->rmesh.density->empty())
  {
    VertexDensity const& vd1(this->rmesh.density[v1]);
    VertexDensity const& vd2(this->rmesh.density[v2]);
    d1 = vd1.valid ? vd1.density : 1.0f;
    d2 = vd2.valid ? vd2.density : 1.0f;
  }

  /* Total amount of density on the edge. */
  float dens = (1.0f / 2.0f) * (d1 + d2) * len;
  float density_left = dens;

  std::size_t verts_inserted = 0;
  float current_len = 0.0f;
  while (max_vertices > 0)
  {
    /* Check if the edge is large enough to handle next teleport. */
    if (density_left < excess_teleport)
    {
        excess_teleport -= density_left;
        return verts_inserted;
    }

    current_len = this->eval_edge_coord(len, d1, d2, current_len, excess_teleport);

#ifdef REMESHER_NAN_CHECKS
    if (std::isnan(current_len))
    {
        // std::cout << "NAN values while resampling, EXIT" << std::endl;
        // std::exit(1);
    }
#endif

    /* Place a sample here. */
    float lambda = current_len / len;
    this->place_feature_sample(v1, v2, lambda);
    verts_inserted += 1;
    max_vertices -= 1;

    /* Consume density and advance. */
    density_left -= excess_teleport;
    excess_teleport = density_spacing;
  }

  return verts_inserted;
}

/* ---------------------------------------------------------------- */

float
ResamplingSampler::eval_edge_coord (float len, float d1, float d2,
    float cur, float dens)
{
  /* Solves with constant x1 for x2 in:
   * \int_x1^x2 d(x) dx == dens
   * where d(x) is the linear density over the edge.
   */

  /* A simpler and more stable version is applied if d1 ~ d2. */
  if (FLOAT_EQ(d1, d2))
  {
    return cur + dens / d1;
  }
  else
  {
    float ddiff = (d2 - d1) / (2.0f * len);
    float x1_const = ddiff * cur * cur + d1 * cur;
    /* ABC-formula here. */
    float x2 = (1.0f / (2.0f * ddiff))
        * (-d1 + std::sqrt(d1 * d1 - 4.0f * ddiff * (-x1_const - dens)));

    return x2;
  }
}

/* ---------------------------------------------------------------- */

void
ResamplingSampler::place_feature_sample (std::size_t v1,
    std::size_t v2, float lambda)
{
  if (EPSILON_EQ(lambda, 0.0f, 0.00001f))
    lambda = 0.0f;
  if (EPSILON_EQ(lambda, 1.0f, 0.00001f))
    lambda = 1.0f;

  if (this->per_edge.empty() || this->per_edge.back().v1 != v1
      || this->per_edge.back().v2 != v2)
    this->per_edge.push_back(SamplesPerEdge(v1, v2, 1));
  else
    this->per_edge.back().samples += 1;
  this->per_edge_pos.push_back(lambda);
}

/* ---------------------------------------------------------------- */

float
ResamplingSampler::get_backbone_density (MeshBackbone const& mb) const
{
  float backbone_density = 0.0f;
  for (BackboneIter j = mb.verts.begin(); j != mb.verts.end(); ++j)
  {
    BackboneIter jp1 = ++BackboneIter(j);
    std::size_t idx1 = *j;
    std::size_t idx2 = (jp1 == mb.verts.end() ? mb.verts.front() : *jp1);

    if (!mb.closed && idx2 == mb.verts.front())
      break;

    backbone_density += this->rmesh.density->get_edge_density(idx1, idx2);
  }
  return backbone_density;
}

/* ================================================================ */

void
ResamplingMesher::exec_meshing (void)
{
  /* Create the new mesh. */
#if RESAMPLING_ON_MESHCOPY
  this->emesh.mesh = this->rmesh.mesh->create_copy();
  this->emesh.features = this->rmesh.features->create_copy();
  this->emesh.vrefs = VertexRefList::create();
  this->emesh.vrefs->init_identity(this->rmesh.mesh, this->rmesh.vinfo);
#else
  this->emesh.mesh = TriangleMesh::create();
  this->emesh.features = FeatureEdges::create();
  this->emesh.vrefs = VertexRefList::create();
#endif

  /* Fetch resampling result. */
  SamplesPerFaceList const& flist = this->sampler->get_samples_per_face();
  SamplesPerEdgeList const& elist = this->sampler->get_samples_per_edge();
  EdgeSamplePositions const& epos = this->sampler->get_edge_sample_pos();

  /* Init delete list for mesh decimation later on. */
  this->dlist.resize(this->emesh.mesh->get_vertices().size(), true);

#if MUTUAL_TESSELLATION_COLORS
  this->emesh.mesh->get_vertex_colors().resize
      (this->emesh.mesh->get_vertices().size(), MUTUAL_COLOR_ORIG);
#endif

  /* Sample the surface and features. */
  this->sample_faces(flist);
  this->sample_features(elist, epos);

#if 0
  // TEST to show mutual tessellation
  this->emesh.features->set_vertex_info
      (VertexInfoList::create(this->emesh.mesh));
  this->emesh.features->set_mesh(this->emesh.mesh);
  return;
#endif


#if RESAMPLING_ON_MESHCOPY
  /* std::cout << "Performing some edge flips..." << std::endl; */
  DelaunayFlips df;
  df.set_mesh(this->emesh.mesh);
  df.set_vertex_info(VertexInfoList::create(this->emesh.mesh));
  df.set_feature_edges(this->emesh.features);
  df.flip_edges();

  if (this->config.perform_decimation)
  {
    // std::cout << "Performing decimation, deleting old vertices..." << std::endl;
    // std::cout << "  Info: The mesh has currently "
    //     << this->emesh.mesh->get_vertices().size() << " vertices" << std::endl;

    MeshDecimation::RelocList reloc;

    MeshDecimation md;
    md.set_mesh(this->emesh.mesh);
    md.set_features(this->emesh.features);
    md.set_delete_list(this->dlist);
    md.set_exact_budget(this->config.sample_amount);
    md.fill_vertex_reloc_list(reloc);
    md.start_decimation();

    /* std::cout << "Decimation finished, fixing vertex references..." << std::endl; */

    /* Fix vertex references. */
    VertexRefListPtr vrefs = VertexRefList::create();
    vrefs->resize(reloc.size());
    for (std::size_t i = 0; i < reloc.size(); ++i)
      vrefs[i] = this->emesh.vrefs[reloc[i]];
    this->emesh.vrefs = vrefs;
  }
  else
  {
    /* Fix feature info. */
    this->emesh.features->set_vertex_info
        (VertexInfoList::create(this->emesh.mesh));
    this->emesh.features->set_mesh(this->emesh.mesh);
  }
#endif
}

/* ---------------------------------------------------------------- */

void
ResamplingMesher::sample_faces (SamplesPerFaceList const& flist)
{
  for (std::size_t i = 0; i < flist.size(); ++i)
    this->sample_triangle(i, flist[i]);
}

/* ---------------------------------------------------------------- */

void
ResamplingMesher::sample_features (SamplesPerEdgeList const& elist,
    EdgeSamplePositions const& epos)
{
#if RESAMPLING_ON_MESHCOPY

  /* Data structures quick access. */
  MeshVertexList& everts = this->emesh.mesh->get_vertices();
  MeshFaceList& efaces = this->emesh.mesh->get_faces();

  /* Create vertex info for feature sampling. */
  VertexInfoListPtr vinfo(VertexInfoList::create(this->emesh.mesh));

  std::size_t pi = 0;
  for (std::size_t i = 0; i < elist.size(); ++i)
  {
    SamplesPerEdge const& spe(elist[i]);

    /* Handle corner vertices. */
    if (spe.v1 == spe.v2)
      this->dlist[spe.v1] = false;

    if (spe.samples == 0)
      continue;

    /* Check if we still operate on valid edges (assertion). */
    if (!vinfo->is_mesh_edge(spe.v1, spe.v2))
    {
      // std::cout << "Error finding edge: " << spe.v1 << "->" << spe.v2
      //     << ", aborting" << std::endl;
      vinfo->debug_vertex(spe.v1);
      vinfo->debug_vertex(spe.v2);
      throw Exception("Error: Expecting edge but cannot find it!");
    }

    /* Find appropirate face to specify vertex references. */
    std::size_t vref_face = MAX_SIZE_T;
    std::size_t vref_edge_off = MAX_SIZE_T;
    bool vref_reverse = false;
    this->rmesh.vinfo->get_face_for_edge(spe.v1, spe.v2,
        &vref_face, &vref_edge_off, &vref_reverse);
    if (vref_face == MAX_SIZE_T)
      throw Exception("Error: Cannot find face for vref");

    /* ---- Create new edge vertices and insert. ---- */

    std::vector<std::size_t> insvec;
    insvec.push_back(spe.v1);

    // Copy vertices because location may change (push_back)
    Vec3f v1(everts[spe.v1]);
    Vec3f v2(everts[spe.v2]);

    for (std::size_t j = 0; j < spe.samples; ++j)
    {
      /* Check if barycentric coordiante was snapped to vertex. */
      if (epos[pi] == 0.0f)
      {
        //std::cout << "Snapped v1! " << spe.v1 << std::endl;
        this->dlist[spe.v1] = false;
#if MUTUAL_TESSELLATION_COLORS
        this->emesh.mesh->get_vertex_colors()[spe.v1] = MUTUAL_COLOR_NEW;
#endif
      }
      else if (epos[pi] == 1.0f)
      {
        //std::cout << "Snapped v2! " << spe.v2 << std::endl;
        this->dlist[spe.v2] = false;
#if MUTUAL_TESSELLATION_COLORS
        this->emesh.mesh->get_vertex_colors()[spe.v2] = MUTUAL_COLOR_NEW;
#endif
      }
      else
      {
        float bary1 = 1.0f - epos[pi];
        float bary2 = epos[pi];
        Vec3f newv(v1 * bary1 + v2 * bary2);
        insvec.push_back(everts.size());
        everts.push_back(newv);
        vinfo->push_back(VertexInfo());

#if MUTUAL_TESSELLATION_COLORS
        this->emesh.mesh->get_vertex_colors().push_back(MUTUAL_COLOR_NEW);
#endif

        /* Provide a vertex refernce. */
        Vec3f bary(0.0f, 0.0f, 0.0f);
        bary[vref_edge_off] = (vref_reverse ? bary2 : bary1);
        bary[(vref_edge_off + 1) % 3] = (vref_reverse ? bary1 : bary2);
        this->emesh.vrefs->push_back(VertexRef
            (vref_face, Vec2f(bary[0], bary[1])));

        /* Update delete list (don't delete that vertex!) */
        this->dlist.push_back(false);
      }

      pi += 1;
    }

    insvec.push_back(spe.v2);

    /* Fixing features. */
    if (!this->emesh.features->empty())
    {
      this->emesh.features->resize(everts.size());
      this->emesh.features->rm_feature_edge(spe.v1, spe.v2);
      for (std::size_t j = 0; j < insvec.size() - 1; ++j)
        this->emesh.features->add_feature_edge(insvec[j], insvec[j+1]);
    }

    /* Insert new triangles every adjacent for edge vertices. */
    VertexInfo::FaceList const& adj_faces = vinfo[spe.v1].adj_faces;
    for (std::size_t j = 0; j < adj_faces.size(); ++j)
    {

      std::size_t face = adj_faces[j];

      /* Find the edge in the triangle, remember direction and 3rd vertex. */
      bool reverse = false;
      std::size_t v3 = MAX_SIZE_T;
      for (std::size_t k = 0; k < 3; ++k)
      {
        std::size_t kp1 = (k + 1) % 3;
        std::size_t adj_v1 = efaces[face * 3 + k];
        std::size_t adj_v2 = efaces[face * 3 + kp1];
        if (spe.v1 == adj_v1 && spe.v2 == adj_v2)
        {
          reverse = false;
          v3 = efaces[face * 3 + (k + 2) % 3];
        }
        else if (spe.v1 == adj_v2 && spe.v2 == adj_v1)
        {
          reverse = true;
          v3 = efaces[face * 3 + (k + 2) % 3];
        }
      }

      if (v3 == MAX_SIZE_T)
        continue;

      /* Now insert new triangles for the adjacent face. */
      //std::cout << "Setting face "  << face << " to "
      //    << insvec[reverse ? 1 : 0] << " " << insvec[reverse ? 0 : 1]
      //    << " " << v3 << std::endl;
      efaces[face * 3 + 0] = (MeshVIndex)insvec[reverse ? 1 : 0];
      efaces[face * 3 + 1] = (MeshVIndex)insvec[reverse ? 0 : 1];
      efaces[face * 3 + 2] = (MeshVIndex)v3;
      for (std::size_t k = 1; k < insvec.size() - 1; ++k)
      {
        /* Fix vertex info for v3. */
        vinfo[v3].adj_faces.push_back(efaces.size() / 3);
        /* Add the new face... */
        efaces.push_back((MeshVIndex)insvec[k + (reverse ? 1 : 0)]);
        efaces.push_back((MeshVIndex)insvec[k + (reverse ? 0 : 1)]);
        efaces.push_back((MeshVIndex)v3);
      }
      vinfo->order_and_classify(v3);

      /* Fix vertex info for v2. */
      vinfo->remove_adjacent_face(spe.v2, face);
      vinfo[spe.v2].adj_faces.push_back(efaces.size() / 3 - 1);
      vinfo->order_and_classify(spe.v2);
    }
  }

#else

  MeshVertexList const& rverts = this->rmesh.mesh->get_vertices();

  std::size_t pi = 0;
  for (std::size_t i = 0; i < elist.size(); ++i)
  {
    SamplesPerEdge const& spe(elist[i]);

    /* Insert edge vertices. */
    for (std::size_t j = 0; j < spe.samples; ++j)
    {
      Vec3f v1(rverts[spe.v1]);
      Vec3f v2(rverts[spe.v2]);
      Vec3f p(v1 * (1.0f - epos[pi]) + v2 * epos[pi]);
      this->emesh.mesh->get_vertices().push_back(p);
#if MUTUAL_TESSELLATION_COLORS
      this->emesh.mesh->get_vertex_colors().push_back(MUTUAL_COLOR_NEW);
#endif
      pi += 1;
    }

    /* Check for corner vertices. */
    if (spe.v1 == spe.v2 && spe.samples == 0)
    {
      Vec3f p(rverts[spe.v1]);
      this->emesh.mesh->get_vertices().push_back(p);
#if MUTUAL_TESSELLATION_COLORS
      this->emesh.mesh->get_vertex_colors().push_back(MUTUAL_COLOR_NEW);
#endif
    }
  }

#endif
}

/* ---------------------------------------------------------------- */
// FIXME: Sample triangles according to density function

void
ResamplingMesher::sample_triangle (std::size_t tri, std::size_t samples)
{
  if (samples == 0)
    return;

  /* 
   * <simlan> I need to randomly sample a triangle using barycentric
   * coordinates. I can do this by generating three random numbers and
   * normalizing them, but this creates a non-uniform distribution for
   * arbitrary triangles. Is there a way to create a uniform distribution?
   *
   * <WimC> simlan: pick uniform random x1, x2 in [0,1] and take the
   * lengths of the intervals between 0, x1, x2, 1 (ordering x1 and x2).
   *
   * <WimC> simlan: let T = {(x1,x2) | 0<=x1,x2<=1 and x1<=x2}} then the
   * map (x1,x2) |-> x1, x2-x1, 1-x2 is affine and onto the triangle
   * {(x,y,z) | 0<=x,y,z<=1 and x+y+z=1}.  same for the transposed triangle
   * x1>=x2. so its a double affine cover.
   */
  MeshVertexList const& rverts = this->rmesh.mesh->get_vertices();
  MeshFaceList const& rfaces = this->rmesh.mesh->get_faces();
  MeshVertexList& everts = this->emesh.mesh->get_vertices();
  MeshFaceList& efaces = this->emesh.mesh->get_faces();

  /* Fetch vertices and density values. */
  std::size_t vidx[3] = { rfaces[tri * 3 + 0],
      rfaces[tri * 3 + 1], rfaces[tri * 3 + 2] };
  Vec3f vert[3] = { rverts[vidx[0]], rverts[vidx[1]], rverts[vidx[2]] };
  float dens[3] = { 1.0f, 1.0f, 1.0f };
  if (!this->rmesh.density->empty())
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      VertexDensity const& vd(this->rmesh.density[vidx[i]]);
      dens[i] = (vd.valid ? vd.density : 1.0f);
    }
  }

  
/*
  std::cout << "  Sampling triangle (" << vidx[0] << "," << vidx[1]
      << "," << vidx[2] << ") with density values "
      << dens[0] << " " << dens[1] << " " << dens[2]
      << " with " << samples << " samples..." << std::endl;
  
*/

#if RESAMPLING_ON_MESHCOPY
  /* Maintain an incremental delaunay triangulation. */
  IncDelaunay incdel(vert[0], vert[1], vert[2]);
  incdel.set_insert_epsilon(0.1f);
#endif

  for (std::size_t i = 0; i < samples; ++i)
  {
#if RESAMPLING_ON_MESHCOPY
    /* Insert the new sample into the incremental delaunay triangulation. */
    bool inserted = false;
    std::size_t insert_cnt = 0;
    do
    {
      //Vec3f bary(this->get_random_bary());
      Vec3f bary(this->get_random_bary(dens[0], dens[1], dens[2]));
      inserted = incdel.insert_sample(Vec2f(bary[0], bary[1]));
      insert_cnt += 1;
    }
    while (!inserted && insert_cnt < 1000);

    if (insert_cnt >= 1000)
      std::cout << "Warning: Insertion failed after 1000 tries" << std::endl;

#else
    //Vec3f bary(this->get_random_bary());
    Vec3f bary(this->get_random_bary(dens[0], dens[1], dens[2]));
    Vec3f p(vert[0] * bary[0] + vert[1] * bary[1] + vert[2] * bary[2]);
    everts.push_back(p);
# if MUTUAL_TESSELLATION_COLORS
      this->emesh.mesh->get_vertex_colors().push_back(MUTUAL_COLOR_NEW);
# endif
#endif
  }

#if RESAMPLING_ON_MESHCOPY

  /*
   * Merge the delaunay triangulation with the evolving mesh. The first
   * face in the Delaunay triangulation replaces the current face in the
   * mesh. Remaining faces are pushed to the end. New vertices are pushed
   * to the end, faces are offsetted with current vertex size.
   */
  std::size_t vsize = everts.size();
  IncDelaunay::IDVertexList const& idvl = incdel.get_vertices();
  HalfEdge::FaceList const& hefl = incdel.get_faces();

  if (hefl.empty())
    throw Exception("Face list is unexpectedly empty!");

  if (idvl.size() < 3)
    throw Exception("Vertex list with less than 3 vertices");

  /* Insert vertices to the mesh. */
  for (std::size_t i = 3; i < idvl.size(); ++i)
  {
    IDVertex const& v(idvl[i]);
    Vec3f newpos = vert[0] * v.bary[0] + vert[1] * v.bary[1]
        + vert[2] * (1.0f - v.bary[0] - v.bary[1]);
    everts.push_back(newpos);
# if MUTUAL_TESSELLATION_COLORS
    this->emesh.mesh->get_vertex_colors().push_back(MUTUAL_COLOR_NEW);
# endif
    this->emesh.vrefs->push_back(VertexRef(tri, Vec2f(v.bary[0], v.bary[1])));
    dlist.push_back(false);
  }

  /* Insert faces to the mesh. */
  for (std::size_t i = 0; i < hefl.size(); ++i)
  {
    std::size_t id1 = hefl[i]->edge->vert->id;
    std::size_t id2 = hefl[i]->edge->next->vert->id;
    std::size_t id3 = hefl[i]->edge->next->next->vert->id;
    if (i == 0)
    {
      efaces[tri * 3 + 0] = (MeshVIndex)(id1 < 3 ? vidx[id1] : id1 + vsize - 3);
      efaces[tri * 3 + 1] = (MeshVIndex)(id2 < 3 ? vidx[id2] : id2 + vsize - 3);
      efaces[tri * 3 + 2] = (MeshVIndex)(id3 < 3 ? vidx[id3] : id3 + vsize - 3);
    }
    else
    {
      efaces.push_back((MeshVIndex)(id1 < 3 ? vidx[id1] : id1 + vsize - 3));
      efaces.push_back((MeshVIndex)(id2 < 3 ? vidx[id2] : id2 + vsize - 3));
      efaces.push_back((MeshVIndex)(id3 < 3 ? vidx[id3] : id3 + vsize - 3));
    }
  }

  /* Update mesh features (new samples was inserted). */
  if (!this->emesh.features->empty())
    this->emesh.features->resize(everts.size());

#endif
}

/* ---------------------------------------------------------------- */

Vec3f
ResamplingMesher::get_random_bary (void) const
{
  /* Create random barycentric coordinate. */
  float x1 = (float)std::rand() / (float)RAND_MAX;
  float x2 = (float)std::rand() / (float)RAND_MAX;
  if (x1 > x2)
    std::swap(x1, x2);
  
  return Vec3f(x1, x2 - x1, 1.0f - x2);
}

/* ---------------------------------------------------------------- */

Vec3f
ResamplingMesher::get_random_bary (float d1, float d2, float d3)
{
  if (d1 == d2 && d1 == d3)
    return this->get_random_bary();

  int numRoots = 0;
  while (numRoots == 0)
  {
    double x1 = (double)std::rand() / (double)RAND_MAX;
    double x2 = (double)std::rand() / (double)RAND_MAX;

	double r1 = 0.0;
	double r2 = 0.0;
	double r3 = 0.0;

	double factor = (d1 + d2 + d3)/6.0;

	// get fist transformed random variate
    if (EPSILON_EQ(d2 + d3, 2.0f * d1, 0.0001f))
    {
//      numRoots = gsl_poly_solve_quadratic(1.0, -2.0, x1, &r1, &r2);
      numRoots = gsl_poly_solve_quadratic(-0.5*d1, d1, -factor*x1, &r1, &r2);
	}
	else
    {
      numRoots = gsl_poly_solve_cubic (3.0 * (d1-d2-d3) / (d2+d3-2*d1),
        3.0 * (d2+d3)/(d2+d3-2*d1), -6.0/(d2+d3-2*d1)*factor*x1, &r1, &r2, &r3);
	}

    if (numRoots == 0)
        continue;

    int correctRoot = 0;
	if ((r1 >= 0.0) && (r1 <=1.0))
		correctRoot += 1;
	if ((numRoots > 1) && (r2 >= 0.0) && (r2 <=1.0))
		correctRoot += 2;
	if ((numRoots > 2) && (r3 >= 0.0) && (r3 <=1.0))
		correctRoot += 4;

    double X = 0.0f;
    switch (correctRoot)
    {
	  case 0: throw Exception("None of the X roots in [0,1], density: "
          + Helpers::get_string(d1, 5) + " " + Helpers::get_string(d2, 5)
          + " " + Helpers::get_string(d3, 5) + " roots: "
          + Helpers::get_string(numRoots) + " factor: " + Helpers::get_string(factor, 5));
      case 1: X = r1; break;
      case 2: X = r2; break;
      case 4: X = r3; break;
      default: throw Exception("Multiple X roots in [0,1]");
    }


    numRoots = 0;
    if (d2 == d3)
    {
      numRoots = 1;
      r1 = 0.5 * ((d2+d3-2*d1)*X*X + 2*(d1-d2-d3) * X + d2 + d3) * x2 / (d1*X - d3*X + d3);
    }
    else
    {
//      numRoots = gsl_poly_solve_quadratic(1.0,
//        2.0 * (d1*X - d3*X + d3) / (d2-d3),
//        -1.0 * ((d2+d3-2*d1)*X*X + 2*(d1-d2-d3)*X + d2 + d3) / (d2-d3) * x2,
//        &r1, &r2);
      numRoots = gsl_poly_solve_quadratic(0.5*(d2-d3),
        d1*X - d3*X + d3,
        -0.5 * ((d2+d3-2*d1)*X*X + 2*(d1-d2-d3)*X + d2 + d3) * x2,
        &r1, &r2);
    }

    if (numRoots == 0)
	    continue;

    correctRoot = 0;
	if ((r1 >= 0.0) && (r1 <= 1.0-X))
       correctRoot += 1;
	if ((numRoots > 1) && (r2 >= 0.0) && (r2 <=1.0-X))
      correctRoot += 2;

    double Y = 0.0f;
    switch (correctRoot)
    {
      case 0: throw Exception("None of the Y roots in [0,1]");
      case 1: Y = r1; break;
      case 2: Y = r2; break;
      default: throw Exception("Multiple Y roots in [0,1]");
    }


    return Vec3f((float)X, (float)Y, 1.0f - (float)X - (float)Y);
  }

  throw Exception("Komisch");
}

/* ---------------------------------------------------------------- */


/* ================================================================ */

void
Resampling::exec_resampling (void)
{
  ::srand(0); // Temp

  this->sampler.set_mesh_data(this->rmesh);
  this->sampler.exec_resampling(this->config.sample_amount);

  this->mesher.set_input_mesh(this->rmesh);
  this->mesher.set_config(this->config);
  this->mesher.set_sampler(this->sampler);
  this->mesher.exec_meshing();

  this->emesh = this->mesher.get_output_mesh();
}

REMESHER_NAMESPACE_END
