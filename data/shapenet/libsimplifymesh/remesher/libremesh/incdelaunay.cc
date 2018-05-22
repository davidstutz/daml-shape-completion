#include <fstream>
#include <iostream>

#include "exception.h"
#include "incdelaunay.h"

REMESHER_NAMESPACE_BEGIN

IncDelaunay::IncDelaunay (Vec3f const& v1, Vec3f const& v2, Vec3f const& v3)
{
  this->ins_eps = 0.0001f;

  /* Define the first two vertices and calc the third one. */
  Vec2f v1x(0.0f, 0.0f);
  Vec2f v2x(1.0f, 0.0f);

  float scale = 1.0f / (v1 - v2).length();
  float len1 = (v1 - v3).length() * scale;
  float len2 = (v2 - v3).length() * scale;
  float v3_x = (1.0f / 2.0f) * (len1 * len1 - len2 * len2 + 1.0f);
  float v3_y = std::sqrt(len1 * len1 - v3_x * v3_x);

  Vec2f v3x(v3_x, v3_y);

  //std::cout << "Input triangle converted to (0,0), (1,0), (" << v3_x
  //   << "," << v3_y << ")" << std::endl;

  /* Insert master triangle. */
  this->verts.push_back(IDVertex(v1x, Vec2f(1.0f, 0.0f)));
  this->verts.push_back(IDVertex(v2x, Vec2f(0.0f, 1.0f)));
  this->verts.push_back(IDVertex(v3x, Vec2f(0.0f, 0.0f)));

  HEFace* face = new HEFace();
  HEVertex* vert0 = new HEVertex(0);
  HEVertex* vert1 = new HEVertex(1);
  HEVertex* vert2 = new HEVertex(2);
  HEEdge* edge2 = new HEEdge(vert0, 0, face, 0);
  HEEdge* edge1 = new HEEdge(vert2, 0, face, edge2);
  HEEdge* edge0 = new HEEdge(vert1, 0, face, edge1);
  edge2->next = edge0;
  face->edge = edge0;
  vert0->edge = edge0;
  vert1->edge = edge1;
  vert2->edge = edge2;

  this->push_back(face);
  this->push_back(vert0);
  this->push_back(vert1);
  this->push_back(vert2);
  this->push_back(edge0);
  this->push_back(edge1);
  this->push_back(edge2);
}

/* ---------------------------------------------------------------- */

bool
IncDelaunay::insert_sample (Vec2f const& bary)
{
  if (this->verts.size() < 3)
    throw Exception("IncDelaunay: Not properly initialized!");

  /* Calc vertex position w.r.t. master triangle. */
  Vec2f pos = this->verts[0].pos * bary[0] + this->verts[1].pos * bary[1]
      + this->verts[2].pos * (1.0f - bary[0] - bary[1]);

  /* Find face that contains "pos". */
  bool inserted = false;
  for (std::size_t i = 0; i < this->he_faces.size(); ++i)
  {
    HEFace* face = this->he_faces[i];
    Vec3f lbary = this->get_bary_for_face(face, pos);
    if (!this->is_inside_face(lbary))
      continue;

    /* Inspect for numeric instabilities. */
    if (EPSILON_EQ(lbary[0], 0.0f, this->ins_eps)
        || EPSILON_EQ(lbary[1], 0.0f, this->ins_eps)
        || EPSILON_EQ(lbary[2], 0.0f, this->ins_eps))
    {
      //std::cout << "Numerically instable sample. Skipping" << std::endl;
      return false;
    }

    this->verts.push_back(IDVertex(pos, bary));
    delaunay_insert(face);
    inserted = true;
    break;
  }

  return inserted;
}


/* ---------------------------------------------------------------- */

Vec3f
IncDelaunay::get_bary_for_face (HEFace const* face, Vec2f const& pos) const
{
  // http://www.blackpawn.com/texts/pointinpoly/default.html
  HEEdge const* he1(face->edge);
  HEEdge const* he2(he1->next);
  HEEdge const* he3(he2->next);

  Vec2f v1(this->verts[he1->vert->id].pos);
  Vec2f v2(this->verts[he2->vert->id].pos);
  Vec2f v3(this->verts[he3->vert->id].pos);

  Vec2f e1(v3 - v1);
  Vec2f e2(v2 - v1);
  Vec2f e3(pos - v1);

  float dot11 = e1 * e1;
  float dot12 = e1 * e2;
  float dot13 = e1 * e3;
  float dot22 = e2 * e2;
  float dot23 = e2 * e3;

  float denom = (dot11 * dot22 - dot12 * dot12);
  float u = (dot22 * dot13 - dot12 * dot23) / denom;
  float v = (dot11 * dot23 - dot12 * dot13) / denom;

  return Vec3f(1.0f - u - v, v, u);
}

/* ---------------------------------------------------------------- */

void
IncDelaunay::delaunay_insert (HEFace* face)
{
  HEEdge* he1(face->edge);
  HEEdge* he2(he1->next);
  HEEdge* he3(he2->next);

  Vec2f v1(this->verts[he1->vert->id].pos);
  Vec2f v2(this->verts[he2->vert->id].pos);
  Vec2f v3(this->verts[he3->vert->id].pos);

  HEVertex* v = new HEVertex(this->verts.size() - 1);
  HEFace* f1 = face;
  HEFace* f2 = new HEFace(he2);
  HEFace* f3 = new HEFace(he3);

  HEEdge* f1e1 = f1->edge;
  HEEdge* f1e3 = new HEEdge(he3->vert, 0, f1, f1e1);
  HEEdge* f1e2 = new HEEdge(v, 0, f1, f1e3);

  HEEdge* f2e1 = f2->edge;
  HEEdge* f2e3 = new HEEdge(he1->vert, f1e2, f2, f2e1);
  HEEdge* f2e2 = new HEEdge(v, 0, f2, f2e3);

  HEEdge* f3e1 = f3->edge;
  HEEdge* f3e3 = new HEEdge(he2->vert, f2e2, f3, f3e1);
  HEEdge* f3e2 = new HEEdge(v, f1e3, f3, f3e3);

  he1->next = f1e2;
  he2->face = f2;
  he2->next = f2e2;
  he3->face = f3;
  he3->next = f3e2;
  f1e3->pair = f3e2;
  f1e2->pair = f2e3;
  f2e2->pair = f3e3;
  v->edge = f1e3;

  this->push_back(v);
  this->push_back(f2);
  this->push_back(f3);
  this->push_back(f1e2);
  this->push_back(f1e3);
  this->push_back(f2e2);
  this->push_back(f2e3);
  this->push_back(f3e2);
  this->push_back(f3e3);

  this->restore_delaunay(he1->pair);
  this->restore_delaunay(he2->pair);
  this->restore_delaunay(he3->pair);
}

/* ---------------------------------------------------------------- */

void
IncDelaunay::restore_delaunay (HEEdge* edge)
{
  if (edge == 0)
    return;

  std::size_t vidx[4];
  vidx[0] = edge->vert->id;
  vidx[1] = edge->next->vert->id;
  vidx[2] = edge->pair->vert->id;
  vidx[3] = edge->pair->next->vert->id;

  float angle;
  {
    Vec2f u(this->verts[vidx[0]].pos - this->verts[vidx[1]].pos);
    Vec2f v(this->verts[vidx[2]].pos - this->verts[vidx[1]].pos);
    float cosangle = (u * v) / std::sqrt(u.qlength() * v.qlength());
    cosangle = std::min(1.0f, std::max(-1.0f, cosangle));
    angle = std::acos(cosangle);
    REMESHER_NAN_CHECK(angle)
  }

  {
    Vec2f u(this->verts[vidx[2]].pos - this->verts[vidx[3]].pos);
    Vec2f v(this->verts[vidx[0]].pos - this->verts[vidx[3]].pos);
    float cosangle = (u * v) / std::sqrt(u.qlength() * v.qlength());
    cosangle = std::min(1.0f, std::max(-1.0f, cosangle));
    angle += std::acos(cosangle);
    REMESHER_NAN_CHECK(angle)
  }

  /* Check if delaunay criterion is violated or not. */
  if (angle <= MY_PI || FLOAT_EQ(angle, MY_PI))
    return;

  //std::cout << "  Flipping edge " << edge_id << " (" << edge.vs << ","
  //    << edge.ve << "," << edge.vl << "," << edge.vr << ")" << std::endl;

  HEEdge* e1(edge->next);
  HEEdge* e2(e1->next);
  HEEdge* pedge(edge->pair);
  HEEdge* e3(pedge->next);
  HEEdge* e4(e3->next);

  /* Update vertices, faces and edges. */
  edge->vert->edge = e1;
  pedge->vert->edge = e3;

  edge->face->edge = edge;
  pedge->face->edge = pedge;

  edge->vert = e3->vert;
  edge->next = e4;
  pedge->vert = e1->vert;
  pedge->next = e2;
  e1->next = edge;
  e2->face = pedge->face;
  e2->next = e3;
  e3->next = pedge;
  e4->face = edge->face;
  e4->next = e1;

  /* Recursively restore delaunay. */
  this->restore_delaunay(e1->pair);
  this->restore_delaunay(e2->pair);
}

/* ---------------------------------------------------------------- */

void
IncDelaunay::write_debug_mesh (std::string const& filename)
{
  std::ofstream out(filename.c_str());
  if (!out.good())
    throw Exception("Error writing debug mesh");

  out << "OFF" << std::endl;
  out << this->verts.size() << " "
      << this->he_faces.size() << " 0" << std::endl;

  for (std::size_t i = 0; i < this->verts.size(); ++i)
  {
    Vec2f const& p(this->verts[i].pos);
    out << p[0] << " " << p[1] << " 0.0" << std::endl;
  }
  for (std::size_t i = 0; i < this->he_faces.size(); ++i)
  {
    HEFace const* f(this->he_faces[i]);
    HEEdge const* e1(f->edge);
    HEEdge const* e2(e1->next);
    HEEdge const* e3(e2->next);
    out << "3 " << e1->vert->id << " " << e2->vert->id << " "
        << e3->vert->id << std::endl;
  }
  out.close();
}

REMESHER_NAMESPACE_END
