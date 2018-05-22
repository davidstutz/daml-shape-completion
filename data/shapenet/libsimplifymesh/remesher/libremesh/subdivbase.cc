#include <iostream>

#include "elapsedtimer.h"
#include "subdivbase.h"

REMESHER_NAMESPACE_BEGIN

void
SubdivBase::start_subdiv (void)
{
  ElapsedTimer t;

  this->subdiv_impl();

  std::cout << "Subdivision took " << t.get_elapsed() << "ms. ";
  std::cout << "New mesh has " << (this->mesh->get_faces().size() / 3)
      << " faces and " << this->mesh->get_vertices().size()
      << " vertices." << std::endl;
}

/* ---------------------------------------------------------------- */

void
SubdivBase::insert_faces (NewVertexIndices const& newvi)
{
  /* Vertex info is invalidated. Clear it to save memory. */
  this->vinfo->clear();

  MeshFaceList& faces = this->mesh->get_faces();
  std::size_t face_amount = faces.size() / 3;

  /* All vertices are inserted now. Triangles are inserted now,
   * subdividing each existing triangle. */
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    std::size_t new_faces = 0;
    for (std::size_t j = 0; j < 3; ++j)
      if (newvi[i * 3 + j] != MAX_SIZE_T)
        new_faces += 1;

    if (new_faces == 0)
      continue;

    if (new_faces == 1)
    {
      /* Find edge with the new vertex. */
      std::size_t eidx = MAX_SIZE_T;
      for (std::size_t j = 0; j < 3; ++j)
        if (newvi[i * 3 + j] != MAX_SIZE_T)
          eidx = j;

      std::size_t v0idx = faces[i * 3 + eidx];
      std::size_t v1idx = newvi[i * 3 + eidx];
      std::size_t v2idx = faces[i * 3 + (eidx + 1) % 3];
      std::size_t v3idx = faces[i * 3 + (eidx + 2) % 3];

      faces[i * 3 + 0] = (MeshVIndex)v0idx;
      faces[i * 3 + 1] = (MeshVIndex)v1idx;
      faces[i * 3 + 2] = (MeshVIndex)v3idx;

      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v2idx);
      faces.push_back((MeshVIndex)v3idx);
    }
    else if (new_faces == 2)
    {
      /* Find edge without the new vertex. */
      std::size_t eidx = MAX_SIZE_T;
      for (std::size_t j = 0; j < 3; ++j)
        if (newvi[i * 3 + j] == MAX_SIZE_T)
          eidx = j;

      std::size_t v0idx = faces[i * 3 + (eidx + 1) % 3];
      std::size_t v1idx = newvi[i * 3 + (eidx + 1) % 3];
      std::size_t v2idx = faces[i * 3 + (eidx + 2) % 3];
      std::size_t v3idx = newvi[i * 3 + (eidx + 2) % 3];
      std::size_t v4idx = faces[i * 3 + eidx];

      faces[i * 3 + 0] = (MeshVIndex)v0idx;
      faces[i * 3 + 1] = (MeshVIndex)v1idx;
      faces[i * 3 + 2] = (MeshVIndex)v4idx;

      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v2idx);
      faces.push_back((MeshVIndex)v3idx);

      faces.push_back((MeshVIndex)v3idx);
      faces.push_back((MeshVIndex)v4idx);
      faces.push_back((MeshVIndex)v1idx);
    }
    else if (new_faces == 3)
    {
      /* Shorthands for the vertex indices. */
      std::size_t v0idx = faces[i * 3 + 0];
      std::size_t v2idx = faces[i * 3 + 1];
      std::size_t v4idx = faces[i * 3 + 2];
      std::size_t v1idx = newvi[i * 3 + 0];
      std::size_t v3idx = newvi[i * 3 + 1];
      std::size_t v5idx = newvi[i * 3 + 2];

      faces[i * 3 + 0] = (MeshVIndex)v0idx;
      faces[i * 3 + 1] = (MeshVIndex)v1idx;
      faces[i * 3 + 2] = (MeshVIndex)v5idx;

      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v2idx);
      faces.push_back((MeshVIndex)v3idx);

      faces.push_back((MeshVIndex)v3idx);
      faces.push_back((MeshVIndex)v4idx);
      faces.push_back((MeshVIndex)v5idx);

      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v3idx);
      faces.push_back((MeshVIndex)v5idx);
    }
  }
}

REMESHER_NAMESPACE_END
