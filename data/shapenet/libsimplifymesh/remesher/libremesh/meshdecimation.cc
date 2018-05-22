#include <cstdlib>
#include <iostream>

#include "exception.h"
#include "averageplane.h"
#include "triangulator.h"
#include "meshcleanup.h"
#include "meshdecimation.h"

REMESHER_NAMESPACE_BEGIN

void
MeshDecimation::start_decimation (void)
{
  /* A few consistency checks... */
  if (this->mesh.get() == 0)
    throw Exception("Mesh not set!");

  if (this->dlist == 0 || this->dlist->empty())
    throw Exception("Delete list not initialized!");

  if (this->dlist->size() != this->mesh->get_vertices().size())
    throw Exception("Invalid delete list given");

  /* Calculate vertex information. */
  this->vinfo = VertexInfoList::create(this->mesh);
  this->vertex_amount = this->vinfo->size();

  /* Try to delete each marked vertex. */
  std::size_t deleted = 0;
  std::vector<std::size_t> retry_list;
  for (std::size_t i = 0; i < this->dlist->size(); ++i)
  {
    if (!this->dlist->at(i))
      continue;

    this->dlist->at(i) = this->delete_vertex(i);
    if (this->dlist->at(i) == false)
      retry_list.push_back(i);
    else
      deleted += 1;
  }

  // std::cout << "Deleted " << deleted << " vertices. retrying..." << std::endl;

  /* Retry to delete. */
  deleted = 0;
  for (std::size_t i = 0; i < retry_list.size(); ++i)
  {
    std::size_t index = retry_list[i];
    this->dlist->at(index) = this->delete_vertex(index);
    if (this->dlist->at(index) == true)
      deleted += 1;
  }
  // std::cout << "Done, deleted " << deleted << " more vertices." << std::endl;

  /* Try to match the specified vertex budget. */
  if (this->vertex_budget != 0)
  {
    while (this->vertex_amount > this->vertex_budget)
    {
      std::size_t rand = std::rand() % this->dlist->size();

      /* Skip already deleted vertices. */
      if (this->dlist->at(rand))
        continue;

      /* Skip feature vertices. */
      if (this->features.get() != 0 && !this->features->empty())
      {
        VertexInfo::VertexList vl;
        this->vinfo->get_adjacent_vertices(rand, vl);
        std::size_t fedges = this->features->count_edges_for_vertex(rand, vl);
        if (fedges > 0)
          continue;
      }

      /* Try to delete vertex. */
      this->dlist->at(rand) = this->delete_vertex(rand);
    }
  }

  /* Clean the mesh from deleted faces and vertices. */
  MeshCleanup cleaner(this->mesh);
  cleaner.set_vertex_info(this->vinfo);
  cleaner.set_vertex_reloc_list(this->reloc);
  cleaner.cleanup_deleted(*this->dlist);

  this->vinfo.reset();

  /* Make a copy of the mesh to free unused space. */
  this->mesh = this->mesh->create_copy();

  /* Update features. */
  this->features->set_mesh(this->mesh);
  this->features->set_vertex_info(VertexInfoList::create(this->mesh));
  this->features->fix_index_relocation(*this->reloc);
}

/* ---------------------------------------------------------------- */

bool
MeshDecimation::delete_vertex (std::size_t index)
{
  /* Check if the vertex should really be deleted. */
  VertexInfo const& vi = this->vinfo[index];
  switch (vi.vclass)
  {
    case VERTEX_CLASS_COMPLEX:
      // std::cout << "Warning: Cannot delete complex vertex!" << std::endl;
      return false;

    case VERTEX_CLASS_UNREFERENCED:
      this->vertex_amount -= 1;
      return true;

    case VERTEX_CLASS_SIMPLE:
      if (vi.adj_faces.size() < 3)
      {
        // std::cout << "Warning: Skipping vertex with <3 faces" << std::endl;
        return false;
      }
      break;

    case VERTEX_CLASS_BORDER:
      if (vi.adj_faces.size() < 2)
        return false;
      break;

    default:
      return false;
  }

  /* Get some information about the vertex. */
  VertexInfo::VertexList vlist;
  this->vinfo->get_adjacent_vertices(index, vlist);

  #if 0
  /* Check for complex vertex loops. */
  if (this->vinfo->is_complex_vertex_loop(vlist))
  {
    std::cout << "Warning: Skipping vertex with complex loop" << std::endl;
    return false;
  }
  #endif

  /* Define a (yet invalid) split line. The resulting hole is forced to
   * be triangulated with the given split line, that is to keep features. */
  SplitLine sline(MAX_SIZE_T, MAX_SIZE_T);

  /* Check if the vertex is on a feature crease. Define a splitline then. */
  bool is_global_feature = false;
  if (this->features.get() != 0 && !this->features->empty())
  {
    FeatureVertexEdges ve = this->features->get_edges_for_vertex(index, vlist);
    if (ve.size() == 2)
    {
      is_global_feature = true;
      for (std::size_t i = 0; i < vlist.size(); ++i)
      {
        if (ve[0] == vlist[i])
          sline.first = i;
        else if (ve[1] == vlist[i])
          sline.second = i;
      }

      if (sline.first == MAX_SIZE_T || sline.second == MAX_SIZE_T)
      {
        sline = SplitLine(MAX_SIZE_T, MAX_SIZE_T);
        // std::cout << "Warning: Splitline was not found!" << std::endl;
        return false;
      }
    }
    else if (ve.size() == 1 || ve.size() > 2)
    {
      return false;
    }
  }

  /* Calculate average plane. */
  Plane3f av_plane;
  {
    MeshVertexList const& verts = this->mesh->get_vertices();
    AveragePlane av;
    for (std::size_t i = 0; i < vlist.size(); ++i)
      av.add_point(verts[vlist[i]]);
    av.add_point(verts[index]);
    av_plane = av.get_average();
  }

  /* Triangulate the vertex loop. */
  MeshFaceList flist;
  Triangulator tri(this->mesh);
  tri.set_vertex_info(this->vinfo);
  bool good = tri.triangulate(vlist, flist, av_plane.n, sline);
  if (!good)
  {
    // std::cout << "Warning: Triangulation failed!" << std::endl;
    return false;
  }

  /* Fix features data structure. */
  if (is_global_feature)
  {
    this->features->rm_feature_edge(index, vlist[sline.first]);
    this->features->rm_feature_edge(index, vlist[sline.second]);
    this->features->add_feature_edge(vlist[sline.first], vlist[sline.second]);
    this->features[index].clear();
  }

  /* Replace faces in the mesh with the new triangulation. */
  VertexInfo::FaceList const& aflist = vi.adj_faces;
  this->replace_faces(aflist, flist);
  this->vertex_amount -= 1;

  return true;
}

/* ---------------------------------------------------------------- */

void
MeshDecimation::replace_faces (VertexInfo::FaceList aflist,
    MeshFaceList const& flist)
{
  if (aflist.size() < flist.size() / 3)
    throw Exception("Cannot insert faces while replacing");

  MeshFaceList& faces = this->mesh->get_faces();

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
    faces[fidx * 3 + offset] = flist[i];
    this->vinfo[flist[i]].adj_faces.push_back(fidx);
  }

  /* Invalidate residual faces by setting vertex indices to 0. */
  for (std::size_t i = flist.size() / 3; i < aflist.size(); ++i)
    for (std::size_t j = 0; j < 3; ++j)
      faces[aflist[i] * 3 + j] = 0;

  /* Update the data structure by reordering the faces. */
  for (std::size_t i = 0; i < flist.size(); ++i)
    this->vinfo->order_and_classify(flist[i]);
}

REMESHER_NAMESPACE_END
