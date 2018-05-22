#include <iostream>
#include <cstdlib>

#include "permute.h"
#include "elapsedtimer.h"
#include "meshoptimize.h"

REMESHER_NAMESPACE_BEGIN

void
MeshOptimize::optimize (std::size_t cache_size)
{
  ElapsedTimer t;

  this->optimize_face_ordering(cache_size);
  this->optimize_vertex_ordering();

  // std::cout << "Triangle mesh optimization took: "
  //     << t.get_elapsed() << "ms." << std::endl;
}

/* ---------------------------------------------------------------- */

void
MeshOptimize::optimize_face_ordering (std::size_t cache_size)
{
  /* We're done if there are no faces. */
  if (this->mesh->get_faces().size() == 0)
    return;

  MeshFaceList& faces = this->mesh->get_faces();
  MeshVertexList& verts = this->mesh->get_vertices();

  /* First, a few data structures are initialized. */
  VertexLiveFaceCounts live;
  VertexCacheTimes cts;
  DeadEndStack dead;
  FaceEmittedFlags emitted;
  MeshFaceList output;

  live.resize(verts.size());
  cts.resize(verts.size(), 0);
  emitted.resize(faces.size() / 3, false);
  output.reserve(faces.size());

  /* Initialize live face counts per vertex. */
  for (std::size_t i = 0; i < verts.size(); ++i)
    live[i] = this->vinfo[i].adj_faces.size();

  /* Pick an arbitrary starting vertex. Confusingly this is called f. */
  std::size_t f = 0;

  /* The cache time stamp. */
  std::size_t ts = cache_size + 1;

  /* A cursor variable. */
  std::size_t cur = 1;

  while (f != MAX_SIZE_T)
  {
    VertexCandidates next;

    /* Process adjacent faces. */
    VertexInfo::FaceList& adj = this->vinfo[f].adj_faces;
    for (std::size_t t = 0; t < adj.size(); ++t)
    {
      /* Skip already emitted faces. */
      if (emitted[adj[t]])
        continue;

      /* For each vertex of the adjacent faces... */
      for (std::size_t v = 0; v < 3; ++v)
      {
        std::size_t vidx = faces[adj[t] * 3 + v];

        output.push_back((MeshVIndex)vidx);
        dead.push(vidx);
        next.insert(vidx);
        live[vidx] -= 1;

        if (ts - cts[vidx] > cache_size)
        {
          cts[vidx] = ts;
          ts = ts + 1;
        }
      }

      emitted[adj[t]] = true;
    }

    f = this->get_next_vertex(cur, cache_size, next, cts, ts, live, dead);
  }

  /* Copy the new face list to the triangle mesh. */
  for (std::size_t i = 0; i < faces.size(); ++i)
    faces[i] = output[i];
}

/* ---------------------------------------------------------------- */

std::size_t
MeshOptimize::get_next_vertex (std::size_t& cur, std::size_t cache_size,
    VertexCandidates const& next, VertexCacheTimes const& cts,
    std::size_t ts, VertexLiveFaceCounts const& live,
    DeadEndStack& dead)
{
  std::size_t n = MAX_SIZE_T;
  std::size_t m = MAX_SIZE_T;
  std::size_t p = MAX_SIZE_T;

  for (VertexCandidates::const_iterator v = next.begin(); v != next.end(); v++)
  {
    if (live[*v] == 0)
      continue;

    p = 0;
    if (ts - cts[*v] + 2 * live[*v] <= cache_size)
      p = ts - cts[*v];

    if (p > m || m == MAX_SIZE_T)
    {
      m = p;
      n = *v;
    }
  }

  if (n == MAX_SIZE_T)
    n = this->skip_dead_end(live, dead, cur);

  return n;
}

/* ---------------------------------------------------------------- */

std::size_t
MeshOptimize::skip_dead_end (VertexLiveFaceCounts const& live,
    DeadEndStack& dead, std::size_t& cur)
{
  MeshVertexList const& verts = this->mesh->get_vertices();

  while (!dead.empty())
  {
    std::size_t dv = dead.top();
    dead.pop();

    if (live[dv] > 0)
      return dv;
  }

  while (cur < verts.size() - 1)
  {
    cur += 1;
    if (live[cur] > 0)
      return cur;
  }

  return MAX_SIZE_T;
}

/* ---------------------------------------------------------------- */

void
MeshOptimize::optimize_vertex_ordering (void)
{
  MeshVertexList& verts = this->mesh->get_vertices();
  MeshFaceList& faces = this->mesh->get_faces();

  VertexIndexRelocList tmp_reloc;
  VertexIndexRelocList* ireloc;

  if (this->irlist != 0)
    ireloc = this->irlist;
  else
    ireloc = &tmp_reloc;

  ireloc->clear();
  ireloc->resize(verts.size(), MAX_SIZE_T);

  std::size_t fill_pos = 0;

  /* Walk over vertices in face order, create index relocation,
   * fix vertex indices of faces and pysically relocate vertex. */
  for (std::size_t i = 0; i < faces.size(); ++i)
  {
    /* Test if the vertex has already been reordered. */
    std::size_t vidx = faces[i];
    if (ireloc->at(vidx) == MAX_SIZE_T)
    {
      /* Vertex is to be reordered. */
      ireloc->at(vidx) = fill_pos;
      fill_pos += 1;
    }

    /* Update the face. */
    faces[i] = (MeshVIndex)ireloc->at(vidx);
  }

  /* Move the remaining unreferenced vertices to the end. */
  for (std::size_t i = 0; fill_pos < ireloc->size() && i < ireloc->size(); ++i)
  {
    if (ireloc->at(i) == MAX_SIZE_T)
    {
      ireloc->at(i) = fill_pos;
      fill_pos += 1;
    }
  }

  /* Permute the vertices using the index permutation. */
  Permute<Vec3f, std::size_t>::permute_reloc(verts, *ireloc);
}

/* ---------------------------------------------------------------- */

void
MeshOptimize::randomize (void)
{
  MeshFaceList& faces = this->mesh->get_faces();

  std::size_t face_amount = faces.size() / 3;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    std::size_t rnd = ::rand() % face_amount;
    std::swap(faces[i * 3 + 0], faces[rnd * 3 + 0]);
    std::swap(faces[i * 3 + 1], faces[rnd * 3 + 1]);
    std::swap(faces[i * 3 + 2], faces[rnd * 3 + 2]);
  }
}

REMESHER_NAMESPACE_END
