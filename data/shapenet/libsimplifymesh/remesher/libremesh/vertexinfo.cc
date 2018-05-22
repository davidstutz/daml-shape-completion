#include <algorithm>
#include <iostream>
#include <set>

#include "exception.h"
#include "vertexinfo.h"

REMESHER_NAMESPACE_BEGIN

void
VertexInfoList::create_data_structure (void)
{
  MeshVertexList const& verts = this->mesh->get_vertices();
  MeshFaceList const& faces = this->mesh->get_faces();
  std::size_t face_amount = faces.size() / 3;

  this->clear();
  this->resize(verts.size());

  /* Create a list of adjacent faces for each vertex. */
  for (std::size_t i = 0; i < face_amount; ++i)
    for (std::size_t j = 0; j < 3; ++j)
      this->at(faces[i * 3 + j]).adj_faces.push_back(i);
}

/* ---------------------------------------------------------------- */

void
VertexInfoList::order_and_classify_all (void)
{
  for (std::size_t i = 0; i < this->size(); ++i)
    this->order_and_classify(i);

  #if 0
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    std::cout << "Classification for " << i << ": "
        << this->at(i).vclass << ", faces: ";
    for (std::size_t j = 0; j < this->at(i).adj_faces.size(); ++j)
      std::cout << this->at(i).adj_faces[j] << " ";
    std::cout << std::endl;
  }
  #endif
}

/* ---------------------------------------------------------------- */

void
VertexInfoList::order_and_classify (std::size_t index)
{
  //std::cout << "Classifying vertex " << index << std::endl;
  MeshFaceList const& faces = this->mesh->get_faces();

  VertexInfo::FaceList adj_faces = this->at(index).adj_faces;
  if (adj_faces.empty())
  {
    /* This is an unreferenced vertex. It needs to be removed. */
    this->at(index).vclass = VERTEX_CLASS_UNREFERENCED;
    return;
  }

  #if 0
  std::cout << "Adjacent faces: ";
  for (std::size_t i = 0; i < adj_faces.size(); ++i)
    std::cout << adj_faces[i] << " ";
  std::cout << std::endl;
  #endif

  VertexInfo::FaceList ordered;
  ordered.push_back(adj_faces.front());
  adj_faces.erase(adj_faces.begin());

  /* Start with a seed triangle and find the
   * starting and the common vertex index. */
  bool forward = true;
  std::size_t fidx = ordered.front() * 3;
  std::size_t starting = MAX_SIZE_T;
  std::size_t common = MAX_SIZE_T;
  for (std::size_t i = 0; i < 3; ++i)
    if (faces[fidx + i] == index)
    {
      starting = faces[fidx + (i + 1) % 3];
      common = faces[fidx + (i + 2) % 3];
    }

  #if 0
  std::cout << "  Starting with face " << ordered.back() << std::endl;
  std::cout << "  Starting / Common: " << starting
       << " / " << common << std::endl;
  #endif

  /* Sanity check, should not happen. */
  if (starting == MAX_SIZE_T || common == MAX_SIZE_T)
    throw Exception("Triangle does not belong!");

  /* Append faces to the front of the ordered list. */
  while (common != starting)
  {
    /* If there are no more adjacent faces and we're not
     * around yet, it's border vertex. */
    if (adj_faces.empty())
    {
      //std::cout << "  Recognised as border vertex!" << std::endl;
      this->at(index).vclass = VERTEX_CLASS_BORDER;
      this->at(index).adj_faces = ordered;
      return;
    }

    /* Walk over all remaining adjacent faces and find new face in order. */
    std::size_t next = MAX_SIZE_T;
    for (std::size_t i = 0; next == MAX_SIZE_T && i < adj_faces.size(); ++i)
    {
      //std::cout << "  Checking face " << adj_faces[i]
      //    << " for " << common << std::endl;
      for (std::size_t j = 0; next == MAX_SIZE_T && j < 3; ++j)
      {
        fidx = adj_faces[i] * 3;
        if (forward && faces[fidx + j] == index
            && faces[fidx + (j + 1) % 3] == common)
        {
          next = faces[fidx + (j + 2) % 3];
          ordered.push_back(adj_faces[i]);
          adj_faces.erase(adj_faces.begin() + i);
        }
        else if (!forward && faces[fidx + j] == index
            && faces[fidx + (j + 2) % 3] == common)
        {
          next = faces[fidx + (j + 1) % 3];
          ordered.insert(ordered.begin(), adj_faces[i]);
          adj_faces.erase(adj_faces.begin() + i);
        }
      }
    }

    /* Check if next common vertex has been found. */
    if (next == MAX_SIZE_T)
    {
      //std::cout << "  Next common vertex not found" << std::endl;
      if (forward)
      {
        /* There is no more triangle to append. Reverse direction. */
        std::size_t temp = starting;
        starting = common;
        common = temp;
        forward = false;
      }
      else
      {
        /* Direction was reversed already. This is non-2-manifold. */
        this->at(index).vclass = VERTEX_CLASS_COMPLEX;
        return;
      }
    }
    else
    {
      common = next;
    }
  }

  /* Adjacent triangles form a loop. If we have no more triangles,
   * it's a simple vertex. Otherwise it's a complex, non-2-manifold. */
  if (adj_faces.empty())
  {
    this->at(index).vclass = VERTEX_CLASS_SIMPLE;
    this->at(index).adj_faces = ordered;
  }
  else
  {
    this->at(index).vclass = VERTEX_CLASS_COMPLEX;
    return;
  }
}

/* ---------------------------------------------------------------- */

void
VertexInfoList::get_adjacent_vertices (std::size_t index,
    VertexInfo::VertexList& list) const
{
  list.clear();

  /* Special handling for complex vertices. */
  if (this->at(index).vclass == VERTEX_CLASS_COMPLEX)
  {
    std::set<std::size_t> verts;

    MeshFaceList const& faces = this->mesh->get_faces();
    VertexInfo::FaceList const& flist = this->at(index).adj_faces;

    for (std::size_t i = 0; i < flist.size(); ++i)
    {
      std::size_t fidx = flist[i] * 3;
      for (int j = 0; j < 3; ++j)
        if (faces[fidx + j] != index)
          verts.insert(faces[fidx + j]);
    }

    /* Add vertices to the list. */
    for (std::set<std::size_t>::iterator iter = verts.begin();
        iter != verts.end(); ++iter)
      list.push_back(*iter);

    return;
  }

  /* Handling for all other vertex types. */
  MeshFaceList const& faces = this->mesh->get_faces();
  VertexInfo::FaceList const& flist = this->at(index).adj_faces;
  for (std::size_t i = 0; i < flist.size(); ++i)
  {
    std::size_t fidx = flist[i] * 3;
    for (int j = 0; j < 3; ++j)
      if (faces[fidx + j] == index)
      {
        list.push_back(faces[fidx + (j + 1) % 3]);
        if (this->at(index).vclass == VERTEX_CLASS_BORDER
            && i == flist.size() - 1)
          list.push_back(faces[fidx + (j + 2) % 3]);
        break;
      }
  }
}

/* ---------------------------------------------------------------- */

bool
VertexInfoList::remove_adjacent_face (std::size_t index, std::size_t face_id)
{
  VertexInfo::FaceList& adj_faces = this->at(index).adj_faces;
  VertexInfo::FaceList::iterator iter;
  iter = std::find(adj_faces.begin(), adj_faces.end(), face_id);
  if (iter != adj_faces.end())
  {
    adj_faces.erase(iter);
    return true;
  }

  return false;
}

/* ---------------------------------------------------------------- */

bool
VertexInfoList::get_edge_info (std::size_t v1, std::size_t v2,
    std::size_t& face1, std::size_t& face2, std::size_t& v0, std::size_t& v3)
{
  face1 = MAX_SIZE_T;
  face2 = MAX_SIZE_T;
  v0 = MAX_SIZE_T;
  v3 = MAX_SIZE_T;

  MeshFaceList const& faces = this->mesh->get_faces();

  /* Walk over adjacent faces for the vertex v1. */
  VertexInfo::FaceList const& fl = this->at(v1).adj_faces;
  std::size_t v1v2_amount = 0;
  std::size_t v2v1_amount = 0;
  for (std::size_t i = 0; i < fl.size(); ++i)
  {
    /* Inspect each face. If the face contains v1 -> v2, it is face1!
     * If the face contains v2 -> v1, it is face2. */
    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t fv1 = faces[fl[i] * 3 + j];
      std::size_t fv2 = faces[fl[i] * 3 + (j + 1) % 3];
      if (fv1 == v1 && fv2 == v2)
      {
        face1 = fl[i];
        v0 = faces[fl[i] * 3 + (j + 2) % 3];
        v1v2_amount += 1;
      }
      if (fv1 == v2 && fv2 == v1)
      {
        face2 = fl[i];
        v3 = faces[fl[i] * 3 + (j + 2) % 3];
        v2v1_amount += 1;
      }
    }
  }

  if (v1v2_amount > 1 || v2v1_amount > 1)
    return false;

  return true;
}

/* ---------------------------------------------------------------- */

void
VertexInfoList::get_face_for_edge (std::size_t v1, std::size_t v2,
    std::size_t* face, std::size_t* edge_off, bool* reverse)
{
  MeshFaceList const& faces = this->mesh->get_faces();
  VertexInfo::FaceList const& adj_faces = this->at(v1).adj_faces;

  for (std::size_t i = 0; i < adj_faces.size(); ++i)
  {
    std::size_t face_id = adj_faces[i];
    for (std::size_t j = 0; j < 3; ++j)
    {
      if (face)
        *face = face_id;
      if (edge_off)
        *edge_off = j;

      std::size_t jp1 = (j + 1) % 3;
      if (faces[face_id * 3 + j] == v1 && faces[face_id * 3 + jp1] == v2)
      {
        if (reverse)
          *reverse = false;
        return;
      }

      if (faces[face_id * 3 + j] == v2 && faces[face_id * 3 + jp1] == v1)
      {
        if (reverse)
          *reverse = true;
        return;
      }
    }
  }

  if (face)
    *face = MAX_SIZE_T;
  if (edge_off)
    *edge_off = MAX_SIZE_T;
}

/* ---------------------------------------------------------------- */

bool
VertexInfoList::is_mesh_edge (std::size_t v1, std::size_t v2) const
{
  MeshFaceList const& faces = this->mesh->get_faces();
  VertexInfo::FaceList const& fl = this->at(v1).adj_faces;

  for (std::size_t i = 0; i < fl.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      std::size_t fv1(faces[fl[i] * 3 + j]);
      std::size_t fv2(faces[fl[i] * 3 + (j + 1) % 3]);

      if ((fv1 == v1 && fv2 == v2) || (fv1 == v2 && fv2 == v1))
        return true;
    }
  }

  return false;
}

/* ---------------------------------------------------------------- */

VertexInfoList::FaceList
VertexInfoList::get_faces_for_edge (std::size_t v1, std::size_t v2)
{
  MeshFaceList const& faces = this->mesh->get_faces();
  VertexInfo::FaceList const& fl = this->at(v1).adj_faces;

  FaceList ret;

  for (std::size_t i = 0; i < fl.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      std::size_t fv1(faces[fl[i] * 3 + j]);
      std::size_t fv2(faces[fl[i] * 3 + (j + 1) % 3]);

      if ((fv1 == v1 && fv2 == v2) || (fv1 == v2 && fv2 == v1))
        ret.push_back(fl[i]);
    }
  }

  return ret;
}

/* ---------------------------------------------------------------- */

bool
VertexInfoList::is_complex_vertex_loop (VertexInfo::VertexList const& list)
{
  if (list.size() < 3)
    throw Exception("Invalid vertex loop with <3 vertices");

  for (std::size_t i = 0; i < list.size() - 2; ++i)
    for (std::size_t j = i + 2; j < list.size(); ++j)
    {
      if (i == 0 && j + 1 == list.size())
        continue;
      if (this->is_mesh_edge(list[i], list[j]))
        return true;
    }

  return false;
}

/* ---------------------------------------------------------------- */

std::size_t
VertexInfoList::get_memory_usage (void) const
{
  /* This is an approximation only. */
  std::size_t s1 = this->capacity() * sizeof(VertexInfo);
  std::size_t s2 = this->size() * sizeof(std::size_t) * 6;
  return s1 + s2;
}

/* ---------------------------------------------------------------- */

void
VertexInfoList::debug_vertex (std::size_t vid) const
{
  std::cout << "Info for vertex " << vid << ", faces: ";
  VertexInfo::FaceList const& fl(this->at(vid).adj_faces);
  MeshFaceList const& faces(this->mesh->get_faces());
  
  for (std::size_t i = 0; i < fl.size(); ++i)
  {
    std::cout << fl[i] << " (";
    for (int j = 0; j < 3; ++j)
        std::cout << faces[fl[i] * 3 + j] << " ";
    std::cout << ") ";
  }
  std::cout << std::endl;
}

REMESHER_NAMESPACE_END
