#include <iostream>

#include "exception.h"
#include "subdivlinear.h"

REMESHER_NAMESPACE_BEGIN

void
SubdivLinear::subdiv_impl (void)
{
  MeshFaceList& faces = this->mesh->get_faces();
  MeshVertexList& verts = this->mesh->get_vertices();

  /* Every edge in the mesh is subdivided and a new vertex in inserted.
   * The indices of the vertices are stored in this temp list. */
  NewVertexIndices newvidx;
  newvidx.resize(faces.size(), MAX_SIZE_T);

  /* Walk over all triangles. Get the three neighboring triangles.
   * Subdivide the three edges if not done yet. Register the
   * new vertices in the temporary list to mark the edge as divided. */
  std::size_t face_amount = faces.size() / 3;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    std::size_t fidx = i * 3;

    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t jp1 = (j + 1) % 3;
      std::size_t v1idx = faces[fidx + j];
      std::size_t v2idx = faces[fidx + jp1];

      std::size_t face1, face2, v0idx, v3idx;
      bool unique = this->vinfo->get_edge_info(v1idx, v2idx,
          face1, face2, v0idx, v3idx);

      if (!unique)
        continue;

      if (face1 == MAX_SIZE_T)
        throw Exception("Subdivision: Face1 was not identified");

      if (face2 != MAX_SIZE_T)
      {
        std::size_t jf2 = MAX_SIZE_T;
        for (std::size_t k = 0; k < 3; ++k)
          if (faces[face2 * 3 + k] == v2idx)
            jf2 = k;

        if (jf2 == MAX_SIZE_T)
          throw Exception("Subdivision: Error finding face2 index");

        if (newvidx[face1 * 3 + j] != newvidx[face2 * 3 + jf2])
          throw Exception("Subdivision: Edge info does not match");

        if (newvidx[face1 * 3 + j] != MAX_SIZE_T)
          continue;

        newvidx[face2 * 3 + jf2] = verts.size();
      }

      newvidx[face1 * 3 + j] = verts.size();
      Vec3f newvert = (verts[v1idx] + verts[v2idx]) / 2.0f;
      verts.push_back(newvert);

      //std::cout << "Inserted new vertex " << verts.size() - 1 << std::endl;
    }
  }

  this->insert_faces(newvidx);
}

REMESHER_NAMESPACE_END
