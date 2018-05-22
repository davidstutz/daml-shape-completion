#include <iostream>

#include "exception.h"
#include "modelloader.h"
#include "meshcleanup.h"
#include "trianglemesh.h"
#include "delaunayflips.h"
#include "meshoptimize.h"
#include "resampling.h"
#include "permute.h"
#include "interface.h"

REMESHER_NAMESPACE_BEGIN

Interface::Interface (void)
{
  this->reference_mesh = TriangleMesh::create();
}

/* ---------------------------------------------------------------- */

void
Interface::load_model (std::string const& modelfile)
{
  this->reference_mesh = ModelLoader::load_model(modelfile);
  this->reference_mesh->ensure_normals();
  this->reference_mesh->memory_debug();

  this->reference_vinfo = VertexInfoList::create();
  this->reference_density = DensityField::create();
  this->reference_features = FeatureEdges::create();

  this->reference_density->set_mesh(this->reference_mesh);
  this->reference_density->set_vertex_info(this->reference_vinfo);
  this->reference_features->set_mesh(this->reference_mesh);
  this->reference_features->set_vertex_info(this->reference_vinfo);
}

/* ---------------------------------------------------------------- */

void
Interface::exec_oversampling (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  /* Release some data that is invalidated anyway. */
  this->reference_mesh->clear_normals();
  this->reference_density->clear();

  /* Perform oversampling. */
  Oversampling os;
  os.set_config(this->oversampling_conf);
  os.set_mesh(this->reference_mesh);
  os.set_vertex_info(this->reference_vinfo);
  os.set_feature_edges(this->reference_features);
  os.start_oversampling();

  this->reference_vinfo->clear();
  this->reference_mesh->ensure_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_feature_extraction (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  this->reference_features->set_config(this->feature_edges_conf);
  this->reference_features->set_mesh(this->reference_mesh);
  this->reference_features->set_vertex_info(this->reference_vinfo);
  this->reference_features->extract_features();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_density_calculation (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  this->reference_density->set_config(this->density_field_conf);
  this->reference_density->set_mesh(this->reference_mesh);
  this->reference_density->set_vertex_info(this->reference_vinfo);
  this->reference_density->set_feature_edges(this->reference_features);
  this->reference_density->calculate_density_field();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_simplification (void)
{
  this->assure_reference_mesh();

  /* Execute the simplification. An index relocation table is filled. */
  Simplification::VertexIndexRelocList reloc;
  Simplification s;
  s.set_config(this->simplification_conf);
  s.set_mesh(this->reference_mesh, true);
  s.set_density_field(this->reference_density);
  s.set_features(this->reference_features, true);
  s.fill_vertex_reloc_list(reloc);
  s.start_simplification();

  /* Take the simplified mesh as evolving mesh to be processed. */
  this->evolving_mesh = s.get_mesh();
  this->evolving_features = s.get_features();

  /* Some information is recalculated for the old and new mesh. */
  this->evolving_mesh->recalc_normals();
  this->evolving_mesh->memory_debug();
  this->evolving_vinfo = VertexInfoList::create(this->evolving_mesh);
  this->reference_vinfo->calc_for_mesh(this->reference_mesh);

  /* Vertex references are initialized using the index relocation. */
  std::size_t new_vertices = this->evolving_mesh->get_vertices().size();
  this->evolving_vref = VertexRefList::create();
  this->evolving_vref->resize(new_vertices);
  for (std::size_t i = 0; i < new_vertices; ++i)
  {
    std::size_t orig_idx = reloc[i];
    VertexInfo::FaceList const& adj_faces
        = this->reference_vinfo[orig_idx].adj_faces;

    if (adj_faces.empty())
      throw Exception("Vertex without adjacent faces");

    MeshFaceList const& faces = this->reference_mesh->get_faces();
    std::size_t face_idx = adj_faces[0];
    std::size_t face_off = MAX_SIZE_T;
    for (int j = 0; j < 3; ++j)
      if (faces[face_idx * 3 + j] == orig_idx)
        face_off = j;

    if (face_off == MAX_SIZE_T)
      throw Exception("Adjacent face without vertex");

    this->evolving_vref[i].face = face_idx;
    this->evolving_vref[i].bary = Vec2f((face_off == 0 ? 1.0f : 0.0f),
        (face_off == 1 ? 1.0f : 0.0f));

    #if 0
    std::cout << "New face: " << face_idx << " vertices: "
      << this->reference_mesh->get_faces()[face_idx * 3 + 0] << " "
      << this->reference_mesh->get_faces()[face_idx * 3 + 1] << " "
      << this->reference_mesh->get_faces()[face_idx * 3 + 2]
      << " Bary: " << this->evolving_vref[i].bary << std::endl;
    #endif
  }

  this->optimize_evolving_mesh();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_lloyd (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();
  this->assure_evolving_mesh();
  this->assure_vertex_references();

  LloydRelaxation l;
  l.set_config(this->lloyd_conf);
  l.set_reference_mesh(this->reference_mesh);
  l.set_reference_info(this->reference_vinfo);
  l.set_evolving_mesh(this->evolving_mesh);
  l.set_evolving_info(this->evolving_vinfo);
  l.set_reference_features(this->reference_features);
  l.set_evolving_features(this->evolving_features);
  l.set_vertex_reflist(this->evolving_vref);
  l.set_reference_density(this->reference_density);
  l.start_relaxation();

  //this->optimize_evolving_mesh();
  this->evolving_mesh->recalc_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_area_equalization (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();
  this->assure_evolving_mesh();
  this->assure_vertex_references();

  AreaBasedRelaxation l;
  l.set_config(this->area_equal_conf);
  l.set_reference_mesh(this->reference_mesh);
  l.set_reference_info(this->reference_vinfo);
  l.set_evolving_mesh(this->evolving_mesh);
  l.set_evolving_info(this->evolving_vinfo);
  l.set_evolving_features(this->evolving_features);
  l.set_reference_features(this->reference_features);
  l.set_vertex_reflist(this->evolving_vref);
  l.set_reference_density(this->reference_density);
  l.start_relaxation();

  this->optimize_evolving_mesh();
  this->evolving_mesh->recalc_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::exec_resampling (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  this->evolving_features.reset();
  this->evolving_vinfo.reset();
  this->evolving_vref.reset();
  this->evolving_mesh.reset();

  /* Perform the resampling. */
  Resampling rs;
  rs.set_config(this->resampling_conf);
  rs.set_reference_mesh(this->reference_mesh);
  rs.set_reference_vinfo(this->reference_vinfo);
  rs.set_reference_features(this->reference_features);
  rs.set_reference_density(this->reference_density);
  rs.exec_resampling();

  /* Aquire the new data. */
  this->evolving_mesh = rs.get_resampled_mesh();
  this->evolving_features = rs.get_resampled_features();
  this->evolving_vref = rs.get_vertex_references();

  /* Recalculate some information. */
  this->evolving_mesh->recalc_normals();
  this->evolving_mesh->memory_debug();
  this->evolving_vinfo = VertexInfoList::create(this->evolving_mesh);
}

/* ---------------------------------------------------------------- */

void
Interface::clean_reference_mesh (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  this->reference_density->clear();
  this->reference_features->clear();

  MeshCleanup mc(this->reference_mesh);
  mc.set_vertex_info(this->reference_vinfo);
  mc.cleanup_mesh(MY_FLT_EPS);

  this->reference_vinfo->clear();
  this->reference_mesh->ensure_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::clean_duplicated_vertices (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  MeshCleanup mc(this->reference_mesh);
  mc.set_vertex_info(this->reference_vinfo);
  mc.cleanup_duplicated(0.0001f);

  this->reference_vinfo->clear();
  this->reference_features->clear();
  this->reference_density->clear();
  this->reference_mesh->recalc_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::copy_reference_mesh (void)
{
  this->assure_reference_mesh();

  this->evolving_mesh = this->reference_mesh->create_copy(false);
  this->evolving_mesh->ensure_normals();
  this->evolving_mesh->memory_debug();

  /* Vertex information is prepared for the old and new mesh. */
  this->evolving_vinfo = VertexInfoList::create(this->evolving_mesh);
  this->reference_vinfo->calc_for_mesh(this->reference_mesh);

  /* Copy the features from the reference mesh. */
  this->evolving_features = this->reference_features->create_copy();
  this->evolving_features->set_mesh(this->evolving_mesh);
  this->evolving_features->set_vertex_info(this->evolving_vinfo);

  /* Initialize the vertex references. */
  this->evolving_vref = VertexRefList::create();
  this->evolving_vref->init_identity(this->evolving_mesh, this->evolving_vinfo);
}

/* ---------------------------------------------------------------- */

void
Interface::copy_evolving_mesh (void)
{
  this->assure_evolving_mesh();

  this->evolving_vref->clear();
  this->reference_vinfo->clear();
  this->reference_density->clear();

  this->reference_mesh.reset();

  this->reference_mesh = this->evolving_mesh->create_copy();
  if (this->evolving_features.get() != 0)
  {
    this->reference_features.reset();
    this->reference_features = this->evolving_features->create_copy();
  }
}

/* ---------------------------------------------------------------- */

void
Interface::optimize_reference_mesh (void)
{
  this->assure_reference_mesh();
  this->assure_reference_vinfo();

  MeshOptimize opt;
  opt.set_mesh(this->reference_mesh);
  opt.set_vertex_info(this->reference_vinfo);
  opt.optimize();

  this->reference_vinfo->clear();
  this->reference_mesh->recalc_normals();
}

/* ---------------------------------------------------------------- */

void
Interface::optimize_evolving_mesh (void)
{
  this->assure_evolving_mesh();

  MeshOptimize::VertexIndexRelocList reloc;
  MeshOptimize opt;
  opt.set_mesh(this->evolving_mesh);
  opt.set_vertex_info(this->evolving_vinfo);
  opt.fill_vertex_reloc_list(reloc);
  opt.optimize();

  this->evolving_mesh->recalc_normals();
  this->evolving_vinfo->calc_for_mesh(this->evolving_mesh);

  /* Now adapt the vertex references. */
  Permute<VertexRef, std::size_t>::permute_reloc(*this->evolving_vref, reloc);

  /* Adapt the feature edges. */
  FeatureEdgesPtr old_features = this->evolving_features->create_copy();
  this->evolving_features->clear();
  this->evolving_features->resize(old_features->size());
  for (std::size_t i = 0; i < old_features->size(); ++i)
    for (std::size_t j = 0; j < old_features[i].size(); ++j)
      this->evolving_features->add_feature_edge
          (reloc[i], reloc[old_features[i][j]]);
}

/* ---------------------------------------------------------------- */

void
Interface::exec_delaunay_flips (void)
{
  this->assure_evolving_mesh();

  if (this->evolving_vinfo->empty())
    this->evolving_vinfo->calc_for_mesh(this->evolving_mesh);

  DelaunayFlips df;
  df.set_max_iterations(25);
  df.set_mesh(this->evolving_mesh);
  df.set_vertex_info(this->evolving_vinfo);
  df.set_feature_edges(this->evolving_features);
  df.flip_edges();
}

/* ---------------------------------------------------------------- */

void
Interface::assure_reference_mesh (void)
{
  if (this->reference_mesh.get() == 0
      || this->reference_mesh->get_vertices().empty())
    throw Exception("Reference mesh is not loaded!");
}

/* ---------------------------------------------------------------- */

void
Interface::assure_reference_vinfo (void)
{
  if (this->reference_vinfo.get() == 0)
    throw Exception("Reference vertex info invalid!");

  if (this->reference_vinfo->empty())
    this->reference_vinfo->calc_for_mesh(this->reference_mesh);
}

/* ---------------------------------------------------------------- */

void
Interface::assure_evolving_mesh (void)
{
  if (this->evolving_mesh.get() == 0
      || this->evolving_mesh->get_vertices().empty())
    throw Exception("Evolving mesh is not loaded");
}

/* ---------------------------------------------------------------- */

void
Interface::assure_vertex_references (void)
{
  if (this->evolving_vref.get() == 0 || this->evolving_vref->size()
      != this->evolving_mesh->get_vertices().size())
    throw Exception("Evolving mesh vertex references not initialized!");
}

/* ---------------------------------------------------------------- */

void
Interface::assure_reference_features (void)
{
  if (this->reference_features.get() == 0 || this->reference_features->size()
      != this->reference_mesh->get_vertices().size())
  {
    this->assure_reference_vinfo();
    this->reference_features = FeatureEdges::create();
    this->reference_features->set_mesh(this->reference_mesh);
    this->reference_features->set_vertex_info(this->reference_vinfo);
    this->reference_features->resize(this->reference_mesh->get_vertices().size());
  }
}

REMESHER_NAMESPACE_END
