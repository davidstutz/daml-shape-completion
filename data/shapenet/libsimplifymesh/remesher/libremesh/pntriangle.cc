#include <iostream>
#include "exception.h"
#include "pntriangle.h"

REMESHER_NAMESPACE_BEGIN

void
PNTriangle::init (std::size_t face_id)
{
  if (this->mesh.get() == 0)
    throw Exception("No mesh set!");

  MeshFaceList const& faces(this->mesh->get_faces());
  MeshVertexList const& verts(this->mesh->get_vertices());
  MeshNormalList const& normals(this->mesh->get_vertex_normals());

  for (int i = 0; i < 3; ++i)
  {
    std::size_t vid(faces[face_id * 3 + i]);
    this->v[i] = verts[vid];
    this->n[i] = normals[vid];
  }

  /* Fix some normals if features are given. */
  if (this->features.get() != 0 && !this->features->empty())
  {
    Vec3f const& fnormal = this->mesh->get_face_normals()[face_id];
    for (int i = 0; i < 3; ++i)
    {
      int ip1 = (i + 1) % 3;
      std::size_t vid_i(faces[face_id * 3 + i]);
      std::size_t vid_ip1(faces[face_id * 3 + ip1]);

      if (this->features->is_feature_edge(vid_i, vid_ip1))
      {
        std::cout << "Setting face normal to vertex normal!" << std::endl;
        this->n[i] = fnormal;
        this->n[ip1] = fnormal;
      }
    }
  }
}

/* ---------------------------------------------------------------- */

Vec3f
PNTriangle::get_position (Vec3f const& bary)
{
  float w12 = (v[1] - v[0]).scalar(n[0]);
  float w21 = (v[0] - v[1]).scalar(n[1]);
  float w23 = (v[2] - v[1]).scalar(n[1]);
  float w32 = (v[1] - v[2]).scalar(n[2]);
  float w31 = (v[0] - v[2]).scalar(n[2]);
  float w13 = (v[2] - v[0]).scalar(n[0]);

  Vec3f const& b300(v[0]);
  Vec3f const& b030(v[1]);
  Vec3f const& b003(v[2]);
  Vec3f b210 = (v[0] * 2.0f + v[1] - n[0] * w12) / 3.0f;
  Vec3f b120 = (v[1] * 2.0f + v[0] - n[1] * w21) / 3.0f;
  Vec3f b021 = (v[1] * 2.0f + v[2] - n[1] * w23) / 3.0f;
  Vec3f b012 = (v[2] * 2.0f + v[1] - n[2] * w32) / 3.0f;
  Vec3f b102 = (v[2] * 2.0f + v[0] - n[2] * w31) / 3.0f;
  Vec3f b201 = (v[0] * 2.0f + v[2] - n[0] * w13) / 3.0f;
  Vec3f E = (b210 + b120 + b021 + b012 + b102 + b201) / 6.0f;
  Vec3f V = (v[0] + v[1] + v[2]) / 3.0f;
  Vec3f b111 = E + (E - V) / 2.0f;

  Vec3f b200 = b300 * bary[0] + b210 * bary[1] + b201 * bary[2];
  Vec3f b110 = b210 * bary[0] + b120 * bary[1] + b111 * bary[2];
  Vec3f b101 = b201 * bary[0] + b111 * bary[1] + b102 * bary[2];
  Vec3f b020 = b120 * bary[0] + b030 * bary[1] + b021 * bary[2];
  Vec3f b011 = b111 * bary[0] + b021 * bary[1] + b012 * bary[2];
  Vec3f b002 = b102 * bary[0] + b012 * bary[1] + b003 * bary[2];

  Vec3f b100 = b200 * bary[0] + b110 * bary[1] + b101 * bary[2];
  Vec3f b010 = b110 * bary[0] + b020 * bary[1] + b011 * bary[2];
  Vec3f b001 = b101 * bary[0] + b011 * bary[1] + b002 * bary[2];

  return b100 * bary[0] + b010 * bary[1] + b001 * bary[2];
}

REMESHER_NAMESPACE_END
