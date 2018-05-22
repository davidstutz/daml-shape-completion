#include <iostream>

#include "elapsedtimer.h"
#include "exception.h"
#include "edgeflips.h"

REMESHER_NAMESPACE_BEGIN

void
EdgeFlipsBase::flip_edges (void)
{
  if (this->mesh.get() == 0)
    throw Exception("EdgeFlips: No mesh has been set");

  if (this->mesh->get_faces().empty())
    return;

  if (this->mesh->get_vertices().size() != this->vinfo->size())
    throw Exception("Vertex list does not match vertex info");

  ElapsedTimer total_t;

  std::size_t iteration = 0;
  std::size_t total_flips = 0;
  this->const_flips_iter = 1;
  this->last_flips = 0;
  while (true)
  {
    ElapsedTimer iter_timer;
    std::size_t flips = this->check_flip_edges();
    total_flips += flips;

    // std::cout << "EdgeFlips: Performed " << flips << " edge flips in "
        // << iter_timer.get_elapsed() << "ms." << std::endl;

    iteration += 1;

    if (flips == 0)
      break;

    if (flips == this->last_flips)
      this->const_flips_iter += 1;
    else
    {
      this->last_flips = flips;
      this->const_flips_iter = 1;
    }

    if (this->max_const_flips_iter != 0
        && this->const_flips_iter >= this->max_const_flips_iter)
      break;

    if (this->max_iter != 0 && iteration >= this->max_iter)
      break;
  }

  // std::cout << "EdgeFlips: Performed " << total_flips
      // << " edge flips in " << iteration << " iterations, in "
      // << total_t.get_elapsed() << "ms." << std::endl;
}

/* ---------------------------------------------------------------- */

std::size_t
EdgeFlipsBase::check_flip_edges (void)
{
  /*
   * Iterate over all faces and all three edges per face. Only edges
   * v1 -> v2 are taken if v1 < v2 (v1 and v2 are vertex indices).
   * An actual edge flip is performed if the minimum angle is increased.
   */

  MeshFaceList const& faces = this->mesh->get_faces();

  std::size_t edge_flips = 0;
  std::size_t face_amount = faces.size() / 3;
  for (std::size_t fi = 0; fi < face_amount; ++fi)
  {
    for (std::size_t i = 0; i < 3; ++i)
    {
      std::size_t ip1 = (i + 1) % 3;
      std::size_t v1idx = faces[fi * 3 + i];
      std::size_t v2idx = faces[fi * 3 + ip1];

      /* Only try to flip for vertex indices v1, v2 with v1 < v2. */
      if (v1idx > v2idx)
        continue;

      /* Only flip if we don't have a feature edge. This leads
       * to a constrained delaunay triangulation. */
      if (!this->features->empty()
          && this->features->is_feature_edge(v1idx, v2idx))
        continue;

      std::size_t face1, face2, v0idx, v3idx;
      bool unique = this->vinfo->get_edge_info(v1idx, v2idx,
          face1, face2, v0idx, v3idx);

      /* If the edge is used by more than two faces, skip it. */
      if (!unique)
        continue;

      /* If we cannot find two adjacent faces, the edge is a border. */
      if (face1 == MAX_SIZE_T || face2 == MAX_SIZE_T)
        continue;

      /* If we found an opposive vertex but it's equal to v0idx, we have
       * two faces with the same vertices (in opposite direction). */
      if (v0idx == v3idx)
        continue;

      /* Try to flip. This returns true if the flip is to be performed. */
      if (this->check_flip_edge(v1idx, v2idx, v0idx, v3idx, face1, face2))
      {
        /* The edge should be flipped. */
        this->edge_flip(v1idx, v2idx, v0idx, v3idx, face1, face2);
        edge_flips += 1;
      }
    }
  }

  return edge_flips;
}

/* ---------------------------------------------------------------- */

void
EdgeFlipsBase::edge_flip (std::size_t v1idx, std::size_t v2idx,
    std::size_t v0idx, std::size_t v3idx,
    std::size_t face1, std::size_t face2)

{
  MeshFaceList& faces = this->mesh->get_faces();

  /* Update vertex info. */
  this->vinfo[v0idx].adj_faces.push_back(face2);
  this->vinfo[v3idx].adj_faces.push_back(face1);

  if (!this->vinfo->remove_adjacent_face(v1idx, face2))
    throw Exception("Cannot find face2 in vertex v1");
  if (!this->vinfo->remove_adjacent_face(v2idx, face1))
    throw Exception("Cannot find face1 in vertex v2");

  /* Perform the edge flip. */
  faces[face1 * 3 + 0] = (MeshVIndex)v0idx;
  faces[face1 * 3 + 1] = (MeshVIndex)v1idx;
  faces[face1 * 3 + 2] = (MeshVIndex)v3idx;
  faces[face2 * 3 + 0] = (MeshVIndex)v0idx;
  faces[face2 * 3 + 1] = (MeshVIndex)v3idx;
  faces[face2 * 3 + 2] = (MeshVIndex)v2idx;

  /* Fix vertex info. */
  this->vinfo->order_and_classify(v0idx);
  this->vinfo->order_and_classify(v1idx);
  this->vinfo->order_and_classify(v2idx);
  this->vinfo->order_and_classify(v3idx);
}

REMESHER_NAMESPACE_END
