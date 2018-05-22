#include <map>
#include <algorithm>
#include <iostream>

#include "exception.h"
#include "matrix3.h"
#include "elapsedtimer.h"
#include "featureedges.h"

REMESHER_NAMESPACE_BEGIN

void
FeatureEdges::extract_features (void)
{
  MeshFaceList const& faces = this->mesh->get_faces();
  std::size_t num_faces = faces.size() / 3;

  if (num_faces != this->mesh->get_face_normals().size())
    throw Exception("Feature extraction: No face normals given!");

  std::size_t edges_detected = 0;

  ElapsedTimer t;

  //this->features = VertexFeatureListPtr(new VertexFeatureList);
  //this->features->resize(this->mesh->get_vertices().size());
  this->clear();
  this->resize(this->mesh->get_vertices().size());

  for (std::size_t i = 0; i < num_faces; ++i)
  {
    std::size_t v1 = faces[i * 3 + 0];
    std::size_t v2 = faces[i * 3 + 1];
    std::size_t v3 = faces[i * 3 + 2];

    edges_detected += (int)this->handle_edge(i, v1, v2);
    edges_detected += (int)this->handle_edge(i, v2, v3);
    edges_detected += (int)this->handle_edge(i, v3, v1);
  }

  std::cout << "Feature extraction took " << t.get_elapsed() << "ms. "
      << edges_detected << " edges identified." << std::endl;
}

/* ---------------------------------------------------------------- */

bool
FeatureEdges::handle_edge (std::size_t face, std::size_t v1, std::size_t v2)
{
  if (v1 == v2)
  {
    std::cout << "Warning: Edge has the same endpoints!" << std::endl;
    return false;
  }

  VertexInfoList::FaceList faces = this->vinfo->get_faces_for_edge(v1, v2);

  /* Check for border edges. */
  if (faces.size() == 1)
  {
    if (this->config.border_edges)
    {
      this->add_feature_edge(v1, v2);
      return true;
    }
    else
      return false;
  }

  /* Check for complex edges. */
  if (faces.size() > 2)
  {
    if (this->config.complex_edges)
    {
      if (this->is_feature_edge(v1, v2))
        return false;
      else
        this->add_feature_edge(v1, v2);
      return true;
    }
    else
      return false;
  }

  if (v1 > v2)
    return false;

  std::size_t face1 = faces[0];
  std::size_t face2 = faces[1];

  if (face != face1 && face != face2)
    throw Exception("Error finding expected face!");

  /* Check for the angle. */
  if (this->config.use_angle)
  {
    MeshNormalList const& fnormals = this->mesh->get_face_normals();
    Vec3f const& f1n = fnormals[face1];
    Vec3f const& f2n = fnormals[face2];

    if (f1n.scalar(f2n) <= this->config.angle)
    {
      this->add_feature_edge(v1, v2);
      return true;
    }
  }

  /* Check for a contour edge. */
  if (this->config.extract_contours)
  {
    MeshNormalList& fn = this->mesh->get_face_normals();
    float at = this->config.angle_theta;
    float ap = this->config.angle_phi;

    Matrix3f rot;
    rot[0] = std::cos(at);
    rot[1] = 0.0f;
    rot[2] = std::sin(at);
    rot[3] = std::sin(at) * std::sin(ap);
    rot[4] = std::cos(ap);
    rot[5] = -std::cos(at) * std::sin(ap);
    rot[6] = -std::sin(at) * std::cos(ap);
    rot[7] = std::sin(ap);
    rot[8] = std::cos(at) * std::cos(ap);

    Vec3f xn = rot * Vec3f(0.0f, 0.0f, -1.0f);
    Vec3f const& n1 = fn[face1];
    Vec3f const& n2 = fn[face2];
    if ((xn * n1 <= 0.0f && xn * n2 > 0.0f)
        || (xn * n1 > 0.0f && xn * n2 <= 0.0f))
    {
      this->add_feature_edge(v1, v2);
      return true;
    }
  }

  return false;
}

/* ---------------------------------------------------------------- */

FeatureVertexEdges
FeatureEdges::get_edges_for_vertex (std::size_t index) const
{
  VertexInfo::VertexList vlist;
  this->vinfo->get_adjacent_vertices(index, vlist);
  return this->get_edges_for_vertex(index, vlist);
}

/* ---------------------------------------------------------------- */

FeatureVertexEdges
FeatureEdges::get_edges_for_vertex (std::size_t index,
    VertexInfo::VertexList const& vlist) const
{
  FeatureVertexEdges ret = this->at(index);
  for (std::size_t i = 0; i < vlist.size(); ++i)
  {
    FeatureVertexEdges const& ve = this->at(vlist[i]);
    for (std::size_t j = 0; j < ve.size(); ++j)
      if (ve[j] == index)
        ret.push_back(vlist[i]);
  }

  return ret;
}

/* ---------------------------------------------------------------- */

std::size_t
FeatureEdges::count_edges_for_vertex (std::size_t index,
    VertexInfo::VertexList const& vlist) const
{
  std::size_t ret = this->at(index).size();
  for (std::size_t i = 0; i < vlist.size(); ++i)
  {
    FeatureVertexEdges const& ve = this->at(vlist[i]);
    for (std::size_t j = 0; j < ve.size(); ++j)
      if (ve[j] == index)
        ret += 1;
  }

  return ret;
}

/* ---------------------------------------------------------------- */

bool
FeatureEdges::is_feature_edge (std::size_t index1, std::size_t index2) const
{
  if (index1 > index2)
    std::swap(index1, index2);

  FeatureVertexEdges const& ve = this->at(index1);

  FeatureVertexEdges::const_iterator iter
      = std::find(ve.begin(), ve.end(), index2);
  return iter != ve.end();
}

/* ---------------------------------------------------------------- */

void
FeatureEdges::add_feature_edge (std::size_t index1, std::size_t index2)
{
  if (index1 > index2)
    std::swap(index1, index2);

  this->at(index1).push_back(index2);
}

/* ---------------------------------------------------------------- */

bool
FeatureEdges::rm_feature_edge (std::size_t index1, std::size_t index2)
{
  if (index1 > index2)
    std::swap(index1, index2);

  FeatureVertexEdges& ve = this->at(index1);
  FeatureVertexEdges::iterator iter = std::find(ve.begin(), ve.end(), index2);
  if (iter != ve.end())
  {
    ve.erase(iter);
    return true;
  }

  return false;
}

/* ---------------------------------------------------------------- */

void
FeatureEdges::fix_index_relocation (FeatureEdges::RelocList const& reloc)
{
  if (this->empty())
    return;

  /* Build a reverse relocation map, mapping old indices to new indices. */
  std::map<std::size_t, std::size_t> rev_reloc;
  for (std::size_t i = 0; i < reloc.size(); ++i)
    rev_reloc[reloc.at(i)] = i;

  FeatureEdgesPtr old_features = this->create_copy();
  this->clear();
  this->resize(reloc.size());
  for (std::size_t i = 0; i < reloc.size(); ++i)
  {
    FeatureVertexEdges const& ve = old_features[reloc.at(i)];
    for (std::size_t j = 0; j < ve.size(); ++j)
    {
      std::size_t new_vidx = rev_reloc[ve[j]];
      this->add_feature_edge(new_vidx, i);
    }
  }
}

/* ---------------------------------------------------------------- */

void
FeatureEdges::debug_dump (void) const
{
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    std::cout << "Vertex " << i << ": ";
    for (std::size_t j = 0; j < this->at(i).size(); ++j)
      std::cout << this->at(i)[j] << ", ";
    std::cout << std::endl;
  }
}

REMESHER_NAMESPACE_END
