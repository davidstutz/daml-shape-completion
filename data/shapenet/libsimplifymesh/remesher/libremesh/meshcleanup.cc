#include <set>
#include <algorithm>
#include <iostream>

#include "elapsedtimer.h"
#include "exception.h"
#include "meshcleanup.h"

REMESHER_NAMESPACE_BEGIN

void
MeshCleanup::cleanup_deleted (VertexDeletedList& dlist)
{
  /* Remove normal information, which is invalidated in any case. */
  this->mesh->clear_normals();

  /* Create some short hands. */
  MeshVertexList& verts = this->mesh->get_vertices();
  MeshFaceList& faces = this->mesh->get_faces();
  MeshVertexColorList& colors = this->mesh->get_vertex_colors();
  bool has_colors = (colors.size() == verts.size());
  std::size_t fsize = faces.size() / 3;

  /* Count the remaining vertices. */
  std::size_t vertex_amount = 0;
  for (std::size_t i = 0; i < dlist.size(); ++i)
    vertex_amount += (int)!dlist[i];

  /* Check if the vertex info is to be created first. */
  if (this->vinfo->empty())
  {
    // std::cout << "Warning: Cacluating vertex info for cleanup!" << std::endl;
    this->vinfo->calc_for_mesh(this->mesh);
  }

  /* If a vertex index relocation list is requested, prepare it. */
  if (this->ireloc != 0)
  {
    this->ireloc->clear();
    this->ireloc->resize(vertex_amount);
    for (std::size_t i = 0; i < vertex_amount; ++i)
      this->ireloc->at(i) = i;
  }

  /* Cleaning of the vertices. */
  std::size_t free_iter = 0;
  std::size_t elem_iter = verts.size();
  while (free_iter < verts.size())
  {
    /* Seach a next free position. */
    while (free_iter < verts.size() && !dlist[free_iter])
      free_iter += 1;

    /* Search next element to move. */
    elem_iter -= 1;
    while (elem_iter > free_iter && dlist[elem_iter])
      elem_iter -= 1;

    if (elem_iter <= free_iter)
      break;

    /* At this point, the element elem_iter is moved to the free element
     * free_iter. Update definition of adjacent faces accordingly. */
    VertexInfo::FaceList& flist = this->vinfo[elem_iter].adj_faces;
    for (std::size_t i = 0; i < flist.size(); ++i)
      for (int j = 0; j < 3; ++j)
        if (faces[flist[i] * 3 + j] == elem_iter)
          faces[flist[i] * 3 + j] = (MeshVIndex)free_iter;

    /* Move the vertex and mark as deleted. */
    verts[free_iter] = verts[elem_iter];
    dlist[elem_iter] = true;
    dlist[free_iter] = false;

    /* If the mesh contains colors, fix these accordingly. */
    if (has_colors)
      colors[free_iter] = colors[elem_iter];

    /* If a vertex index reloc list is requested, log that index relocation. */
    if (this->ireloc != 0)
      this->ireloc->at(free_iter) = elem_iter;
  }

  /* Finally shrink the vertex vector. */
  if (free_iter < verts.size())
  {
    verts.resize(free_iter);
    if (has_colors)
      colors.resize(free_iter);
  }

  /* We don't need additional information anymore. */
  dlist.clear();
  this->vinfo->clear();

#if 0
  /* Sanity check to check if invalid vertex indices have been created. */
  for (std::size_t i = 0; i < fsize; ++i)
    for (std::size_t j = 0; j < 3; ++j)
      if (faces[i * 3 + j] >= verts.size())
        std::cout << "Warning: Invalid vertex ID " << faces[i * 3 + j]
            << " detected!" << std::endl;
  return;
#endif

  /* Cleaning of the faces. */
  free_iter = 0;
  elem_iter = fsize;
  while (free_iter < fsize)
  {
    /* Search the next deleted face. */
    while (free_iter < fsize
        && (faces[free_iter * 3 + 0] != 0
        || faces[free_iter * 3 + 1] != 0
        || faces[free_iter * 3 + 2] != 0))
      free_iter += 1;

    /* Search the next valid face. */
    elem_iter -= 1;
    while (elem_iter > free_iter
        && faces[elem_iter * 3 + 0] == 0
        && faces[elem_iter * 3 + 1] == 0
        && faces[elem_iter * 3 + 2] == 0)
      elem_iter -= 1;

    if (elem_iter <= free_iter)
      break;

    /* Move the face. */
    for (int i = 0; i < 3; ++i)
    {
      faces[free_iter * 3 + i] = faces[elem_iter * 3 + i];
      faces[elem_iter * 3 + i] = 0;
    }
  }

  if (free_iter < fsize)
    faces.resize(free_iter * 3);

  if (vertex_amount != verts.size())
    throw Exception("Cleanup: Invalid vertex count detected!");
}

/* ---------------------------------------------------------------- */

void
MeshCleanup::cleanup_mesh (float epsilon)
{
  MeshVertexList const& verts = this->mesh->get_vertices();
  if (verts.empty())
    return;

  ElapsedTimer timer;

  /* Check if the vertex info is to be created first. */
  if (this->vinfo->empty())
  {
    // std::cout << "Warning: Cacluating vertex info for cleanup!" << std::endl;
    this->vinfo->calc_for_mesh(this->mesh);
  }

  /* Create a vector to track deleted vertices. */
  VertexDeletedList dlist;
  dlist.resize(verts.size(), false);

  /* Merge vertices that are very close. */
  std::size_t deleted_vertices = 0;
  std::size_t deleted_faces = 0;
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    /* Skip vertex if it is already deleted. */
    if (dlist[i])
      continue;

    if (this->vinfo[i].vclass == VERTEX_CLASS_UNREFERENCED)
    {
      // std::cout << "Deleting unreferenced vertex " << i << std::endl;
      dlist[i] = true;
      deleted_vertices += 1;
      continue;
    }

    if (this->vinfo[i].vclass == VERTEX_CLASS_SIMPLE
        && this->vinfo[i].adj_faces.size() < 3)
    {
      // std::cout << "Warning: Simple vertex " << i << " with "
          // << this->vinfo[i].adj_faces.size()
          // << " < 3 adjacent faces!" << std::endl;
    }

    VertexInfo::VertexList vlist;
    this->vinfo->get_adjacent_vertices(i, vlist);

    for (std::size_t j = 0; j < vlist.size(); ++j)
    {
      if (i == vlist[j])
      {
        /* Vertex is its own neighbour. */
        std::size_t invalid_face = this->vinfo[i].adj_faces[j];
        // std::cout << "Deleting invalid face "
            // << invalid_face << "... " << std::flush;
        bool deleted = this->delete_invalid_face(invalid_face);
        // std::cout << (deleted ? "success" : "failed") << std::endl;
        deleted_faces += (int)deleted;
      }
      else if (!dlist[vlist[j]]
          && EPSILON_EQ(verts[i][0], verts[vlist[j]][0], epsilon)
          && EPSILON_EQ(verts[i][1], verts[vlist[j]][1], epsilon)
          && EPSILON_EQ(verts[i][2], verts[vlist[j]][2], epsilon))
      {
        /* Vertex very near to an adjacent vertex. */
        // std::cout << "Merging vertices " << i << " and " << vlist[j]
            // << "... " << std::flush;
        std::size_t deleted = this->merge_vertices(i, vlist[j]);
        dlist[vlist[j]] = (deleted > 0);
        deleted_vertices += (int)(deleted > 0);
        deleted_faces += deleted;
        // std::cout << (deleted ? "success" : "failed") << std::endl;
      }
    }
  }

  /* Clean deleted vertices and faces. */
  this->cleanup_deleted(dlist);

  // std::cout << "Cleaning mesh took " << timer.get_elapsed() << "ms. "
      // << "Deleted " << deleted_vertices << " vertices "
      // << "and " << deleted_faces << " faces." << std::endl;
}

/* ---------------------------------------------------------------- */

bool
MeshCleanup::delete_invalid_face (std::size_t face)
{
  MeshFaceList& faces = this->mesh->get_faces();

  /* Reset face and collect vertices to be fixed. */
  std::set<std::size_t> verts;
  for (int i = 0; i < 3; ++i)
  {
    verts.insert(faces[face * 3 + i]);
    faces[face * 3 + i] = 0;
  }

  for (std::set<std::size_t>::iterator iter = verts.begin();
      iter != verts.end(); iter++)
  {
    while (this->vinfo->remove_adjacent_face(*iter, face)) { }
  }

  return true;
}

/* ---------------------------------------------------------------- */

std::size_t
MeshCleanup::merge_vertices (std::size_t a, std::size_t b)
{
  MeshFaceList& faces = this->mesh->get_faces();
  MeshVertexList& verts = this->mesh->get_vertices();
  VertexInfo::FaceList& adj_a = this->vinfo[a].adj_faces;
  VertexInfo::FaceList& adj_b = this->vinfo[b].adj_faces;

  /* Get common triangles (set intersection of adj_a and adj_b). */
  std::set<std::size_t> common_tris;
  for (std::size_t i = 0; i < adj_a.size(); ++i)
    for (std::size_t j = 0; j < adj_b.size(); ++j)
      if (adj_a[i] == adj_b[j])
        common_tris.insert(adj_a[i]);

  if (common_tris.empty())
    return 0;

  /* Fix triangles adjacent to b (to be adjacent to a). */
  for (std::size_t i = 0; i < adj_b.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
      if (faces[adj_b[i] * 3 + j] == b)
        faces[adj_b[i] * 3 + j] = (MeshVIndex)a;
  }

  /* Delete common triangles and update the third vertex. */
  for (std::set<std::size_t>::iterator iter = common_tris.begin();
      iter != common_tris.end(); iter++)
  {
    while (this->vinfo->remove_adjacent_face(a, *iter)) { }
    while (this->vinfo->remove_adjacent_face(b, *iter)) { }

    for (int j = 0; j < 3; ++j)
    {
      std::size_t fidx = *iter * 3 + j;
      if (faces[fidx] != a && faces[fidx] != b)
      {
        while (this->vinfo->remove_adjacent_face(faces[fidx], *iter)) {}
        this->vinfo->order_and_classify(faces[fidx]);
      }
      faces[fidx] = 0;
    }
  }

  /* Add remaining faces of b to the face list of a. */
  adj_a.insert(adj_a.end(), adj_b.begin(), adj_b.end());
  adj_b.clear();
  this->vinfo->order_and_classify(a);

  /* Reposition a to be in the center of a and b. */
  verts[a] = (verts[a] + verts[b]) / 2.0f;

  return common_tris.size();
}

/* ---------------------------------------------------------------- */

void
MeshCleanup::cleanup_duplicated (float epsilon)
{
  ElapsedTimer timer;

  Remesher::MeshVertexList& verts = this->mesh->get_vertices();
  Remesher::MeshFaceList& faces = this->mesh->get_faces();

  std::vector<std::size_t> reloc;
  reloc.resize(verts.size());
  for (std::size_t i = 0; i < verts.size(); ++i)
    reloc[i] = i;

  std::size_t relocations = 0;
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    // if (i % 100 == 0)
    //   std::cout << "Iteration " << i << " of " << verts.size()
    //       << ", working..." << std::endl;

    if (reloc[i] != i)
      continue;

    for (std::size_t j = i + 1; j < verts.size(); ++j)
    {
      if (reloc[j] != j)
        continue;

      Remesher::Vec3f& v1 = verts[i];
      Remesher::Vec3f& v2 = verts[j];
      if (EPSILON_EQ(v1[0], v2[0], epsilon)
          && EPSILON_EQ(v1[1], v2[1], epsilon)
          && EPSILON_EQ(v1[2], v2[2], epsilon))
      {
        reloc[j] = i;
        relocations += 1;
      }
    }
  }

  // std::cout << "Relocating " << relocations
  //     << " vertex indices and fixing faces..." << std::endl;

  std::size_t face_amount = faces.size() / 3;
  std::size_t invalidated = 0;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    MeshVIndex& v1idx = faces[i * 3 + 0];
    MeshVIndex& v2idx = faces[i * 3 + 1];
    MeshVIndex& v3idx = faces[i * 3 + 2];

    v1idx = (MeshVIndex)reloc[v1idx];
    v2idx = (MeshVIndex)reloc[v2idx];
    v3idx = (MeshVIndex)reloc[v3idx];

    /* Invalidate face if it has at least two identical vertices. */
    if (v1idx == v2idx || v1idx == v3idx || v2idx == v3idx)
    {
      v1idx = 0;
      v2idx = 0;
      v3idx = 0;
      invalidated += 1;
    }
  }

  // if (invalidated > 0)
  //   std::cout << "Relocation invalidated " << invalidated << " faces." << std::endl;

  // std::cout << "Deleting invalidated faces and relocated vertices..." << std::endl;

  VertexDeletedList dlist;
  dlist.resize(verts.size());
  for (std::size_t i = 0; i < verts.size(); ++i)
    dlist[i] = (reloc[i] != i);

  this->vinfo->calc_for_mesh(this->mesh);
  this->cleanup_deleted(dlist);

  // std::cout << "Done removing duplicated vertices, took "
  //     << timer.get_elapsed() << " ms." << std::endl;
}

REMESHER_NAMESPACE_END
