#include <iostream>
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

#define MESH_AWPN_NORMALS 1

void
TriangleMesh::recalc_normals (bool face, bool vertex)
{
  if (!face && !vertex)
    return;

  if (face)
  {
    this->face_normals.clear();
    this->face_normals.reserve(this->faces.size() / 3);
  }

  if (vertex)
  {
    this->vertex_normals.clear();
    this->vertex_normals.resize(this->vertices.size(), Vec3f(0.0f, 0.0f, 0.0f));
  }

#if MESH_AWPN_NORMALS

  /*
   * Calculate angle weighted pseudo normals by weighted
   * averaging adjacent face normals.
   */
  for (std::size_t i = 0; i < this->faces.size(); i += 3)
  {
    std::size_t ia = this->faces[i + 0];
    std::size_t ib = this->faces[i + 1];
    std::size_t ic = this->faces[i + 2];

    /* Vertices to calculate the face normal. */
    Vec3f& a = this->vertices[ia];
    Vec3f& b = this->vertices[ib];
    Vec3f& c = this->vertices[ic];
    Vec3f fn = (b - a).cross(c - a).norm();
    REMESHER_NAN_CHECK(fn[0]);

    if (face)
      this->face_normals.push_back(fn);

    if (vertex)
    {
      /* Edge vectors to calculate the angle. */
      Vec3f ab = vertices[ib] - vertices[ia];
      Vec3f bc = vertices[ic] - vertices[ib];
      Vec3f ca = vertices[ia] - vertices[ic];

      /* Caching the length. */
      float abl = ab.length();
      float bcl = bc.length();
      float cal = ca.length();

      /* ...and calculate! */
      this->vertex_normals[ia] += fn * acosf(ab.scalar(-ca) / (abl * cal));
      this->vertex_normals[ib] += fn * acosf((-ab).scalar(bc) / (abl * bcl));
      this->vertex_normals[ic] += fn * acosf(ca.scalar(-bc) / (cal * bcl));
    }
  }

#else

  /*
   * Calculate simple vertex normals by averaging the adjacent face normals.
   */
  for (std::size_t i = 0; i < this->faces.size(); i += 3)
  {
    Vec3f& a = this->vertices[this->faces[i + 0]];
    Vec3f& b = this->vertices[this->faces[i + 1]];
    Vec3f& c = this->vertices[this->faces[i + 2]];
    Vec3f fn = (b - a).cross(c - a);

    if (face)
      this->face_normals.push_back(fn.norm());

    if (vertex)
    {
      this->vertex_normals[this->faces[i + 0]] += fn;
      this->vertex_normals[this->faces[i + 1]] += fn;
      this->vertex_normals[this->faces[i + 2]] += fn;
    }
  }

#endif

  if (vertex)
    for (std::size_t i = 0; i < this->vertex_normals.size(); ++i)
      if (this->vertex_normals[i].qlength() != 0.0f)
        this->vertex_normals[i].norm_self();
}

/* ---------------------------------------------------------------- */

void
TriangleMesh::ensure_normals (bool face, bool vertex)
{
  bool need_face_recalc = (face
      && this->face_normals.size() != this->faces.size() / 3);
  bool need_vertex_recalc = (vertex
      && this->vertex_normals.size() != this->vertices.size());
  this->recalc_normals(need_face_recalc, need_vertex_recalc);
}

/* ---------------------------------------------------------------- */

void
TriangleMesh::invert_faces (void)
{
  std::size_t face_num = this->faces.size() / 3;
  for (std::size_t i = 0; i < face_num; ++i)
  {
    std::size_t idx = i * 3;
    std::size_t tmp = this->faces[idx + 1];
    this->faces[idx + 1] = (MeshVIndex)this->faces[idx + 2];
    this->faces[idx + 2] = (MeshVIndex)tmp;
  }
  this->recalc_normals(true, true);
}

/* ---------------------------------------------------------------- */

void
TriangleMesh::scale_and_center (void)
{
  if (this->vertices.empty())
    return;

  /* Prepare to translate model to the center. */
  Vec3f min(this->vertices[0]);
  Vec3f max(this->vertices[0]);
  for (std::size_t i = 0; i < this->vertices.size(); ++i)
    for (std::size_t j = 0; j < 3; ++j)
    {
      if (min[j] > this->vertices[i][j])
        min[j] = this->vertices[i][j];
      if (max[j] < this->vertices[i][j])
        max[j] = this->vertices[i][j];
    }

  Vec3f move((min + max) / 2.0);

  /* Prepare scaling of model to fit in the einheits cube. */
  float dx = max[0] - min[0];
  float dy = max[1] - min[1];
  float dz = max[2] - min[2];

  float max_yz = MY_MAX(dy, dz);
  float max_xyz = MY_MAX(dx, max_yz);

  /* Translate and scale. */
  for (std::size_t i = 0; i < this->vertices.size(); ++i)
    this->vertices[i] = (this->vertices[i] - move) / max_xyz;
}

/* ---------------------------------------------------------------- */

void
TriangleMesh::memory_debug (void) const
{
  std::size_t vsize = this->vertices.size() * sizeof(Vec3f);
  std::size_t fsize = this->faces.size() * sizeof(MeshVIndex);
  std::size_t fnsize = this->face_normals.size() * sizeof(Vec3f);
  std::size_t vnsize = this->vertex_normals.size() * sizeof(Vec3f);

  // std::cout << "Model memory: V: " << vsize << " F: " << fsize
      // << " VN: " << vnsize << " FN: " << fnsize
      // << " Total: " << (vsize + fsize + fnsize + vnsize) / 1024 << " KB"
      // << std::endl;
}

/* ---------------------------------------------------------------- */

std::size_t
TriangleMesh::get_memory_usage (void) const
{
  // TODO Use capacity
  std::size_t s_verts = this->vertices.capacity() * sizeof(Vec3f);
  std::size_t s_faces = this->faces.capacity() * sizeof(MeshVIndex);
  std::size_t s_vnorm = this->vertex_normals.capacity() * sizeof(Vec3f);
  std::size_t s_fnorm = this->face_normals.capacity() * sizeof(Vec3f);

  return s_verts + s_faces + s_vnorm + s_fnorm;
}

REMESHER_NAMESPACE_END
