#include <limits>
#include <algorithm>
#include <iostream>

#include "exception.h"
#include "elapsedtimer.h"
#include "straightline3.h"
#include "plane.h"
#include "averageplane.h"
#include "meshcleanup.h"
#include "triangulator.h"
#include "simplification.h"

REMESHER_NAMESPACE_BEGIN

Simplification::Simplification (void)
{
  this->ireloc = 0;
}

/* ---------------------------------------------------------------- */

void
Simplification::set_mesh (TriangleMeshPtr mesh, bool copy_mesh)
{
  this->vinfo.reset();
  this->dinfo.clear();

  /* Create a copy of the mesh. Copy vertices, faces and face normals. */
  if (copy_mesh)
  {
    std::cout << "Copying triangle mesh for simplification..." << std::endl;
    this->mesh = mesh->create_copy(false);
    if (this->config.check_local_features)
      this->mesh->ensure_normals(true, false);
  }
  else
  {
    this->mesh = mesh;
  }
}

/* ---------------------------------------------------------------- */

void
Simplification::set_features (FeatureEdgesPtr features, bool copy_features)
{
  if (copy_features)
  {
    this->features = features->create_copy();
  }
  else
  {
    this->features = features;
  }
}

/* ---------------------------------------------------------------- */

void
Simplification::start_simplification (void)
{
  /* Cache amount of vertices and check finish criterion. */
  this->vertex_amount = this->mesh->get_vertices().size();
  if (this->vertex_amount <= this->config.vertex_limit)
  {
    /* Execute this to get a relocation list. (FIXME dirty) */
    this->prepare_vertex_info();
    this->clean_mesh();
    return;
  }

  ElapsedTimer t;

  /* Create a data structure from vertices to faces. */
  this->prepare_vertex_info();

  /* Run the actual simplification algorithm. */
  this->run_simplification();

  /* Clean data structure from deleted faces and vertices. */
  this->clean_mesh();

  std::cout << "Simplification took " << t.get_elapsed() << "ms" << std::endl;
}

/* ---------------------------------------------------------------- */

void
Simplification::prepare_vertex_info (void)
{
  /* The VertexInfoList class creates the basic data structure.
   * The initial classification for each vertex is valid during
   * the whole algorithm, since simplifications will not change
   * any vertex classification. */
  this->vinfo = VertexInfoList::create(this->mesh);

  /* Info if a vertex has been deleted, defaults to false. */
  this->dinfo.resize(this->vinfo->size(), false);
  this->memory_debug();

  /* Point the feature edges to the mesh that is about to be modified. */
  this->features->set_mesh(this->mesh);
  this->features->set_vertex_info(this->vinfo);
}

/* ---------------------------------------------------------------- */

void
Simplification::run_simplification (void)
{
  /* The algorithm works in serveral passes and iterates over all vertices
   * during a pass. After each pass the decimation criterion is loosened.
   */
  float threshold = this->config.initial_threshold;

  /* Each pass may indicate an interuption, that means that there
   * is no more candidate for removal. */
  bool interupted = false;
  std::size_t iteration = 0;
  while (!interupted
      && this->vertex_amount > this->config.vertex_limit
      && threshold < this->config.maximum_threshold)
  {
    std::cout << "Starting pass " << iteration << " with threshold "
        << threshold << "..." << std::endl;
    interupted = this->run_pass(threshold);
    threshold *= this->config.threshold_factor;
    iteration += 1;
  }

  std::cout << "Decimation " << (interupted ? "interupted" : "finished")
      << " after " << iteration << " iterations." << std::endl;
}

/* ---------------------------------------------------------------- */

bool
Simplification::run_pass (float threshold)
{
  //std::cout << "Staring pass with plane-thres " << dist_plane_thres
  //    << " and edge-thres " << dist_edge_thres << std::endl;

  std::size_t num_candidates = 0;
  std::size_t num_removals = 0;

  for (std::size_t i = 0; i < this->vinfo->size(); ++i)
  {
    /* Of course, skip already deleted vertices. */
    if (this->dinfo[i])
      continue;

    switch (this->vinfo[i].vclass)
    {
      /* Always skip complex vertices. */
      default:
      case VERTEX_CLASS_COMPLEX:
        continue;

      /* Skip unreferenced vertices and mark them for removal. */
      case VERTEX_CLASS_UNREFERENCED:
        this->dinfo[i] = true;
        this->vertex_amount -= 1;
        num_removals += 1;
        continue;

      /* Simple, border and inner edge vertices are all handled now. */
      case VERTEX_CLASS_SIMPLE:
      case VERTEX_CLASS_BORDER:
        if (this->handle_removal_candidate(i, threshold))
        {
          this->vertex_amount -= 1;
          num_removals += 1;
        }
        break;
    };

    /* Count this vertex as candidate since we did not skip it. */
    num_candidates += 1;

    /* Check if vertex limit is reached. */
    if (this->vertex_amount <= this->config.vertex_limit)
      break;
  }

  std::cout << "Pass had " << num_candidates << " candidates and removed "
      << num_removals << " vertices, "
      << this->vertex_amount << " vertices left." << std::endl;

  return num_candidates == 0;
}

/* ---------------------------------------------------------------- */

bool
Simplification::handle_removal_candidate (std::size_t index, float thres)
{
  /* Check some preconditions:
   * - Simple vertices must have 3 or more adjacent faces
   * - Border vertices must have 2 or more adjacent faces
   * - Inner edges are treated as simple for this check
   * - Complex and corner vertices are not allowed.
   */
  VertexInfo const& vi = this->vinfo[index];
  switch (vi.vclass)
  {
    case VERTEX_CLASS_SIMPLE:
      if (vi.adj_faces.size() < 3)
      {
        std::cout << "Warning: Skipping vertex with <3 faces." << std::endl;
        return false;
      }
      break;

    case VERTEX_CLASS_BORDER:
      if (vi.adj_faces.size() < 2)
      {
        //std::cout << "Warning: Skipping border with <2 faces." << std::endl;
        return false;
      }
      break;

    default:
      //std::cout << "Warning: Skipping special vertex." << std::endl;
      return false;
  }

  /* Retrieve the list of adjacent vertices. */
  VertexInfo::VertexList vlist;
  this->vinfo->get_adjacent_vertices(index, vlist);

  MeshVertexList const& verts = this->mesh->get_vertices();

  /* Calculate the average plane. The plane is needed for:
   * - Vertex-to-plane distance checks (simple vertices)
   * - Creation of split planes for the triangulation
   */
  AveragePlane av;
  for (std::size_t i = 0; i < vlist.size(); ++i)
    av.add_point(verts[vlist[i]]);
  av.add_point(verts[index]);
  Plane3f av_plane = av.get_average();

  /* Area checks are done here. */
  if (this->config.perform_area_checks)
  {
    if (!this->check_candidate_area(index, thres))
      return false;
  }

  /* Define a splitline, which is invalid for now. Global features
   * might define a splitline, also the fidelity check might define
   * one if local feature checks are enabled. */
  SplitLine sline(0, 0);

  /* Check for global features. This might supply a valid splitline,
   * forbid the removal or just do nothing but allow further processing. */
  if (!this->check_global_features(index, sline, vlist))
    return false;

  bool is_global_feature = false;
  if (sline.first != sline.second)
    is_global_feature = true;

  /* Fidelity checks are done here. The fidelity check makes suggestions
   * for the split line if it detects local inner edge features. */
  if (this->config.perform_fidelity_checks)
  {
    if (!this->check_candidate_fidelity(index, thres, vlist, av_plane, sline))
      return false;
  }

  #if 0
  /* Check if the vertex loop around the vertex to be removed
   * is not complex, i.e. it does not contain edges between loop
   * vertices that are not adjacent. Check is expensive.
   * The check is only executed if there is no splitline defined.
   */
  if (this->vinfo->is_complex_vertex_loop(vlist))
  {
    std::cout << "Warning: Complex vertex loop found." << std::endl;
    return false;
  }
  #endif

  /* Create a new triangulation with the adjacent vertices. */
  MeshFaceList flist;
  Triangulator tri(this->mesh);
  tri.set_vertex_info(this->vinfo);
  bool good = tri.triangulate(vlist, flist, av_plane.n, sline);
  if (!good)
  {
    /* Triangulation failed. */
    //std::cout << "Warning: Triangulation failed" << std::endl;
    return false;
  }

  /* The vertex will be removed. If global feature checks lead to
   * this splitline, fix the feature information. */
  if (is_global_feature)
  {
    this->features->rm_feature_edge(index, vlist[sline.first]);
    this->features->rm_feature_edge(index, vlist[sline.second]);
    this->features->add_feature_edge(vlist[sline.first], vlist[sline.second]);
    this->features[index].clear();
  }

  /* Replace list old faces with the new ones. */
  VertexInfo::FaceList const& aflist = this->vinfo[index].adj_faces;
  this->replace_faces(aflist, flist);

  /* Mark vertex as deleted. */
  this->dinfo[index] = true;

  return true;
}

/* ---------------------------------------------------------------- */

bool
Simplification::check_global_features (std::size_t index, SplitLine& sline,
    VertexInfo::VertexList const& vlist)
{
  /* If there are no global features defined, allow deletion. */
  if (this->features->empty())
    return true;

  /* Check for global features using the feature edges. Corner
   * vertices are never removed. Also inner-edge vertices with
   * neighboring feature edges are not allowed and not removed. */
  FeatureVertexEdges ve = this->features->get_edges_for_vertex(index, vlist);

  /* Forbid removal if global features preservation is requested. */
  if (this->config.keep_global_features && !ve.empty())
    return false;

  if (ve.size() >= 3 || ve.size() == 1)
  {
    /* This is a corner vertex. Disallow simplification. */
    return false;
  }
  else if (ve.size() == 2)
  {
    /* If our vertex is a border vertex, we use regular border
     * simiplification if we have two feature edges and these
     * edges belong to the border. */
    //if (this->vinfo[index].vclass == VERTEX_CLASS_BORDER
    //    && ((vlist.front() == ve[0] && vlist.back() == ve[1])
    //    || (vlist.front() == ve[1] && vlist.back() == ve[0])))
    //  return true;

    /* Create the splitline. The splitline is expressed with
     * adjacent indices, but we have vertex indices. Convert it. */
    SplitLine splitline(MAX_SIZE_T, MAX_SIZE_T);
    for (std::size_t i = 0; i < vlist.size(); ++i)
      if (ve[0] == vlist[i])
        splitline.first = i;
      else if (ve[1] == vlist[i])
        splitline.second = i;

    if (splitline.first == MAX_SIZE_T || splitline.second == MAX_SIZE_T)
    {
      /*
      std::cout << "Warning: Splitline was not found for vertex "
          << index << std::endl;
      std::cout << "  Vertex loop: ";
      for (std::size_t k = 0; k < vlist.size(); ++k)
        std::cout << vlist[k] << " ";
      std::cout << ", VE: " << ve[0] << " " << ve[1] << std::endl;
      std::cout << "  Features for vertex " << index << ": ";
      for (std::size_t k = 0; k < this->features[index].size(); ++k)
        std::cout << this->features[index][k] << " ";
      std::cout << std::endl;
      */

      return false;
    }

    sline = splitline;
    return true;
  }

  /* If there are no global feature edges for that vertex, allow deletion. */
  return true;
}

/* ---------------------------------------------------------------- */

bool
Simplification::check_candidate_fidelity (std::size_t index, float thres,
    VertexInfo::VertexList const& vlist, Plane3f const& plane, SplitLine& sline)
{
  VertexInfo const& vi = this->vinfo[index];
  MeshVertexList const& verts = this->mesh->get_vertices();
  Vec3f const& v = verts[index];

  /* Checks for local features. If it finds index to be a local
   * inner edge vertex, it is treated as a global inner edge vertex. */
  if (this->config.check_local_features
      && sline.first == sline.second
      && vi.vclass == VERTEX_CLASS_SIMPLE)
  {
    std::size_t collected = 0;
    SplitLine splitline;
    VertexInfo::FaceList const& fl = vi.adj_faces;
    for (std::size_t i = 0; i < fl.size(); ++i)
    {
      std::size_t f1idx = fl[i];
      std::size_t f2idx = fl[(i + 1) % fl.size()];
      Vec3f const& f1n = this->mesh->get_face_normals()[f1idx];
      Vec3f const& f2n = this->mesh->get_face_normals()[f2idx];
      if (f1n.scalar(f2n) <= this->config.local_feature_angle)
      {
        switch (collected)
        {
          case 0: splitline.first = (i + 1) % fl.size(); break;
          case 1: splitline.second = (i + 1) % fl.size(); break;
          default: break;
        }
        collected += 1;
      }
    }
    
    //std::cout << "Checking local feature for vertex "
    //    << index << ", collected: " << collected << std::endl;

    /* Vertex has been detected to be locally inner edge.
     * Neighboring vertices are not allowed. */
    if (collected == 2
        && (splitline.first + 1) % vlist.size() != splitline.second
        && (splitline.second + 1) % vlist.size() != splitline.first)
    {
      /* Found a good splitline. Store it. */
      sline = splitline;
    }
  }

  switch (vi.vclass)
  {
    /* A typical simple vertex is checked using distance-to-plane.
     * But it can be simplified as inner edge vertex if a
     * local feature has been detected. */
    case VERTEX_CLASS_SIMPLE:
    {
      if (sline.first != sline.second)
      {
        /* Check distance-to-feature-edge for inner edge vertices. */
        StraightLine3f line(verts[vlist[sline.first]],
            verts[vlist[sline.second]]);
        float dist = line.point_dist(v);
        if (dist > thres * this->config.fidelity_check_penalty)
          return false;
      }
      else
      {
        /* Check distance-to-plane for simple vertices. */
        float dist = std::fabs(plane.point_dist(v));
        if (dist > thres * this->config.fidelity_check_penalty)
          return false;
      }
    }
    break;

    /* Check distance-to-edge for border vertices. */
    case VERTEX_CLASS_BORDER:
    {
      StraightLine3f line(verts[vlist.front()], verts[vlist.back()]);
      float dist = line.point_dist(v);
      if (dist > thres * this->config.fidelity_check_penalty)
        return false;
    }
    break;

    /* We refuse the check the fidelity of unusual vertices. */
    default:
      return false;
  }

  return true;
}

/* ---------------------------------------------------------------- */

bool
Simplification::check_candidate_area (std::size_t index, float thres)
{
  MeshVertexList const& verts = this->mesh->get_vertices();
  MeshFaceList const& faces = this->mesh->get_faces();

  /* Walk over adjacent triangles and sum the areas. */
  float area = 0;
  VertexInfo::FaceList const& flist = this->vinfo[index].adj_faces;
  for (std::size_t i = 0; i < flist.size(); ++i)
  {
    std::size_t a_idx = faces[flist[i] * 3 + 0];
    std::size_t b_idx = faces[flist[i] * 3 + 1];
    std::size_t c_idx = faces[flist[i] * 3 + 2];

    Vec3f const& a = verts[a_idx];
    Vec3f const& b = verts[b_idx];
    Vec3f const& c = verts[c_idx];

    float l1 = (b - a).length();
    float l2 = (c - b).length();
    float l3 = (a - c).length();
    float peri2 = (l1 + l2 + l3) / 2.0f;
    float qtarea = peri2 * (peri2 - l1) * (peri2 - l2) * (peri2 - l3);
    qtarea = MY_MAX(0.0f, qtarea);
    float tarea = std::sqrt(qtarea);

    /* If there is no density function specified on the mesh, we apply a
     * weighing of 1/3 for each vertex wich results in the ordinary area
     * calculation for triangles. If there is a density function specified,
     * each vertex gets a weight of one third of the density at this
     * vertex. */
    float density_factor = 1.0f;
    if (!this->vdens->empty() && this->vdens[index].valid)
    {
      density_factor = 0.0f;
      density_factor += (this->vdens[a_idx].valid
          ? this->vdens[a_idx].density : this->vdens[index].density);
      density_factor += (this->vdens[b_idx].valid
          ? this->vdens[b_idx].density : this->vdens[index].density);
      density_factor += (this->vdens[c_idx].valid
          ? this->vdens[c_idx].density : this->vdens[index].density);
      density_factor /= 3.0f;
      //density_factor = std::pow(density_factor + 0.5f, 5.0f);
    }

    /* The triangle density integrated over the triangle is
     * approximated by multiplying the triangle area with the density
     * factor. The "vertex area" is calculated by dividing the integrated
     * density with 3 to get the area corresponding to one vertex only. */
    area += tarea * density_factor / 3.0f;
  }

  if (area > thres * this->config.area_check_penalty)
    return false;

  return true;
}

/* ---------------------------------------------------------------- */

void
Simplification::replace_faces (VertexInfo::FaceList aflist,
    MeshFaceList const& flist)
{
  if (aflist.size() < flist.size() / 3)
    throw Exception("Cannot insert faces while replacing");

  MeshFaceList& faces = this->mesh->get_faces();

  #if 0
  std::cout << "Adjacent faces: ";
  for (std::size_t i = 0; i < aflist.size(); ++i)
    std::cout << aflist[i] << " ";
  std::cout << std::endl;
  #endif

  /* Remove all references in the vertex info to the old faces. */
  for (std::size_t i = 0; i < aflist.size(); ++i)
  {
    std::size_t fidx = aflist[i];
    for (std::size_t j = 0; j < 3; ++j)
      this->vinfo->remove_adjacent_face(faces[fidx * 3 + j], fidx);
  }

  /* Insert references in the vertex info to the new faces. */
  for (std::size_t i = 0; i < flist.size(); ++i)
  {
    std::size_t offset = i % 3;
    std::size_t fidx = aflist[i / 3];

    #if 0
    if (i % 3 == 0)
    {
      std::cout << "INSERTING: Writing vertex to face " << fidx << ": "
          << flist[i] << " " << flist[i+1] << " " << flist[i+2] << std::endl;
    }
    #endif

    faces[fidx * 3 + offset] = flist[i];

    /* Fix the face normal if local feature checks are performed. */
    if (this->config.check_local_features && i % 3 == 2)
    {
      Vec3f const& a = this->mesh->get_vertices()[flist[i - 2]];
      Vec3f const& b = this->mesh->get_vertices()[flist[i - 1]];
      Vec3f const& c = this->mesh->get_vertices()[flist[i - 0]];
      this->mesh->get_face_normals()[fidx] = (b - a).cross(c - a).norm();
    }

    //std::cout << "  Inserting face " << fidx
    //    << " to vertex " << flist[i] << std::endl;
    this->vinfo[flist[i]].adj_faces.push_back(fidx);
  }

  /* Invalidate residual faces by setting vertex indices to 0. */
  for (std::size_t i = flist.size() / 3; i < aflist.size(); ++i)
    for (std::size_t j = 0; j < 3; ++j)
      faces[aflist[i] * 3 + j] = 0;

  /* Check if new vertex information is properly defined.
   * Update the data structure by reordering the faces. */
  for (std::size_t i = 0; i < flist.size(); ++i)
  {
    VertexClass old_class = this->vinfo[flist[i]].vclass;
    this->vinfo->order_and_classify(flist[i]);
    VertexClass new_class = this->vinfo[flist[i]].vclass;

    std::size_t adj_faces_cnt = this->vinfo[flist[i]].adj_faces.size();
    if (adj_faces_cnt < 3 && new_class != VERTEX_CLASS_BORDER
        && new_class != VERTEX_CLASS_COMPLEX)
    {
      std::cout << "Warning: Vertex with " << adj_faces_cnt
          << " < 3 faces created ["
          << old_class << " -> " << new_class << "]" << std::endl;
    }
  }
}

/* ---------------------------------------------------------------- */

void
Simplification::clean_mesh (void)
{
  std::cout << "Cleaning deleted vertices and faces..." << std::endl;

  /* Deleted vertices and faces are removed from the data structure. */
  MeshCleanup cleaner(this->mesh);
  cleaner.set_vertex_info(this->vinfo);
  cleaner.set_vertex_reloc_list(this->ireloc);
  cleaner.cleanup_deleted(this->dinfo);

  /* Make a copy of the mesh to free unused space. This is because
   * the cleanup algorithm only uses resize on the containers. */
  this->mesh = this->mesh->create_copy();

  /* The set of feature edges for each vertex need to be fixed. */
  this->features->fix_index_relocation(*this->ireloc);
  #if 0
  if (!this->features->empty())
  {
    /* Build a reverse relocation map, mapping old indices to new indices. */
    std::map<std::size_t, std::size_t> rev_ireloc;
    for (std::size_t i = 0; i < this->ireloc->size(); ++i)
      rev_ireloc[this->ireloc->at(i)] = i;

    FeatureEdgesPtr old_features = this->features->create_copy();
    this->features->clear();
    this->features->resize(this->ireloc->size());
    for (std::size_t i = 0; i < this->ireloc->size(); ++i)
    {
      FeatureVertexEdges const& ve = old_features[this->ireloc->at(i)];
      for (std::size_t j = 0; j < ve.size(); ++j)
      {
        std::size_t new_vidx = rev_ireloc[ve[j]];
        this->features->add_feature_edge(new_vidx, i);
      }
    }
  }
  #endif

  std::cout << "Model cleaned. " << this->mesh->get_vertices().size()
      << " vertices and " << this->mesh->get_faces().size() / 3
      << " faces." << std::endl;
}

/* ---------------------------------------------------------------- */

void
Simplification::memory_debug (void)
{
  std::size_t dsize = this->mesh->get_vertices().size() / 8;
  std::size_t isize = sizeof(VertexInfo) * this->mesh->get_vertices().size();
  std::size_t fsize = sizeof(std::size_t) * this->mesh->get_faces().size();

  std::cout << "Simplification memory: V-Info: " << (dsize + isize) / 1000
      << " KB, V-Faces: " << fsize / 1000
      << " KB, Total: " << (isize + dsize + fsize) / 1000
      << " KB" << std::endl;
}

REMESHER_NAMESPACE_END
