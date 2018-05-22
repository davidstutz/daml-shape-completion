#include <iostream>

#include "exception.h"
#include "vertexref.h"

REMESHER_NAMESPACE_BEGIN

void
VertexRefList::init_identity (TriangleMeshPtr mesh, VertexInfoListPtr vinfo)
{
  if (mesh->get_vertices().size() != vinfo->size())
    throw Exception("Vertex info does not match with the mesh!");

  this->clear();
  this->resize(vinfo->size());

  for (std::size_t i = 0; i < vinfo->size(); ++i)
  {
    VertexInfo::FaceList const& adj_faces = vinfo[i].adj_faces;
    if (adj_faces.empty())
    {
      std::cout << "Warning: No info for unreferenced vertex" << std::endl;
      continue;
    }

    std::size_t face = adj_faces[0];
    std::size_t off = MAX_SIZE_T;
    for (std::size_t j = 0; j < 3; ++j)
      if (mesh->get_faces()[face * 3 + j] == i)
        off = j;

    if (off == MAX_SIZE_T)
      throw Exception("Error: Adjacent face without vertex");

    this->at(i).face = face;
    this->at(i).bary = Vec2f(off == 0 ? 1.0f : 0.0f, off == 1 ? 1.0f : 0.0f);
  }
}

REMESHER_NAMESPACE_END
