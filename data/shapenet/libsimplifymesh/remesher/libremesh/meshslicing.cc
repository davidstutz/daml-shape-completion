#include <iostream>
#include <vector>

#include "vec3.h"
#include "exception.h"
#include "triangulator.h"
#include "meshslicing.h"

REMESHER_NAMESPACE_BEGIN

#define SNAP_DIST 0.005f
#define SNAP_ANGLE 0.2f /* 90° - acos(value) = degree eps. */

void
MeshSlicing::lazy_slice (Plane3f const& plane)
{
  MeshVertexList const& verts = this->mesh->get_vertices();
  MeshFaceList const& faces = this->mesh->get_faces();

  std::vector<float> dist;
  dist.resize(verts.size());

  /* Calculate distance for each point to the plane. */
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    dist[i] = plane.point_dist(verts[i]);
    //std::cout << "Distances: " << verts[i] << " " << dist[i] << std::endl;
  }

  /* Walk over all triangles and check check for intersections. */
  std::size_t faces_size = faces.size() / 3;
  for (std::size_t i = 0; i < faces_size; ++i)
  {
    std::size_t fidx = i * 3;
    /* Walk over all edges of the triangle. */
    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t jp1 = (j + 1) % 3;
      std::size_t v1idx = faces[fidx + j];
      std::size_t v2idx = faces[fidx + jp1];
      float dist1 = dist[v1idx];
      float dist2 = dist[v2idx];

      /* Continue if both points are on the same side of the plane. */
      if ((dist1 < 0.0f && dist2 < 0.0f) || (dist1 > 0.0f && dist2 > 0.0f))
        continue;

      /* Continue if one of the points snaps to the plane. */
      if (EPSILON_EQ(dist1, 0.0f, SNAP_DIST)
          || EPSILON_EQ(dist2, 0.0f, SNAP_DIST))
        continue;

      /* Check angle between edge and plane normal. */
      float cosangle = plane.n * (verts[v2idx] - verts[v1idx]).norm();
      if (EPSILON_EQ(cosangle, 0.0f, SNAP_ANGLE))
        continue;

      /* Inserts the new vertex to the edge if not already present. */
      this->process_intersection(v1idx, v2idx, plane);
    }
  }
}

/* ---------------------------------------------------------------- */

void
MeshSlicing::process_intersection (std::size_t v1idx, std::size_t v2idx,
    Plane3f const& plane)
{
  if (v1idx == v2idx)
    throw Exception("Error calculating intersection for the same vertices");

  if (v1idx > v2idx)
    std::swap(v1idx, v2idx);

  MeshVertexList& verts = this->mesh->get_vertices();

  /* Calculate intersection. */
  Vec3f is;
  {
    Vec3f const& v1 = verts[v1idx];
    Vec3f d = verts[v2idx] - v1;
    float lambda = (-v1 * plane.n + plane.d) / (d * plane.n);
    is = v1 + d * lambda;
  }

  /* Check if a similar vertex is already present. */
  EdgeVertices& ev = this->edge_info[std::make_pair(v1idx, v2idx)];
  for (std::size_t i = 0; i < ev.size(); ++i)
  {
    Vec3f const& v = verts[ev[i]];
    if (EPSILON_EQ(is[0], v[0], SNAP_DIST)
        && EPSILON_EQ(is[1], v[1], SNAP_DIST)
        && EPSILON_EQ(is[2], v[2], SNAP_DIST))
      return;
  }

  /* A similar vertex is not present. Insert the vertex. */
  std::size_t new_idx = verts.size();
  verts.push_back(is);
  ev.push_back(new_idx);
}

/* ---------------------------------------------------------------- */

/* A comparison operator that orders vertices on the edge. */
struct EdgeVertexCmp
{
  MeshVertexList const* verts;
  std::size_t vidx;

  EdgeVertexCmp (MeshVertexList const& verts, std::size_t vidx)
    : verts(&verts), vidx(vidx)
  {
  }

  bool operator() (std::size_t v1, std::size_t v2)
  {
    return ((verts->at(vidx) - verts->at(v1)).qlength()
        < (verts->at(vidx) - verts->at(v2)).qlength());
  }
};

/* Triangulation of a triangle using the edge info. */
void
MeshSlicing::triangulate (void)
{
  MeshVertexList const& verts = this->mesh->get_vertices();
  MeshFaceList& faces = this->mesh->get_faces();

  /* Fist, we sort the vertices on all edges in the mesh. */
  for (EdgeMap::iterator iter = this->edge_info.begin();
      iter != this->edge_info.end(); ++iter)
  {
    EdgeVertices& ev = iter->second;
    EdgeVertexCmp cmp(verts, iter->first.first);
    std::sort(ev.begin(), ev.end(), cmp);
  }

  /* Each triangle is triangulated now using the sorted vertices. */
  Triangulator tri(this->mesh);
  std::size_t face_amount = faces.size() / 3;
  std::size_t face_amount100 = face_amount / 100 + 1;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    if (i % face_amount100 == 0)
      std::cout << "Processing face " << i << " ("
          << (int)(100.0f * (float)i / (float)face_amount)
          << "%)" << std::endl;

    std::size_t fidx = i * 3;
    std::size_t v1idx = faces[fidx + 0];
    std::size_t v2idx = faces[fidx + 1];
    std::size_t v3idx = faces[fidx + 2];
    Vec3f const& v1 = verts[v1idx];
    Vec3f const& v2 = verts[v2idx];
    Vec3f const& v3 = verts[v3idx];

    Vec3f ap_normal = (v2 - v1).cross(v3 - v1);

    /* Create a polygon with indices that is going to be triangulated. */
    VertexInfo::VertexList vlist;
    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t a_idx = faces[fidx + j];
      std::size_t b_idx = faces[fidx + (j + 1) % 3];

      vlist.push_back(a_idx);

      bool swapped = false;
      if (a_idx > b_idx)
      {
        std::swap(a_idx, b_idx);
        swapped = true;
      }

      EdgeVertices& ev = this->edge_info[std::make_pair(a_idx, b_idx)];
      for (std::size_t k = 0; k < ev.size(); ++k)
      {
        std::size_t idx = (swapped ? ev.size() - 1 - k : k);
        vlist.push_back(ev[idx]);
      }
    }

    /*
    std::cout << "Vertex loop: ";
    for (std::size_t j = 0; j < vlist.size(); ++j)
      std::cout << vlist[j] << " ";
    std::cout << std::endl;
    */

    /* Triangulate the polygon. */
    MeshFaceList flist;
    Triangulator::SplitLine splitline(0, 0);
    bool good = tri.triangulate(vlist, flist, ap_normal, splitline);
    if (!good)
    {
      std::cout << "Error triangulating!" << std::endl;
      continue;
    }

    /* Update the face list. */
    if (flist.size() < 3)
    {
      std::cout << "Error, no faces in face list" << std::endl;
      continue;
    }

    for (std::size_t j = 0; j < 3; ++j)
      faces[fidx + j] = flist[j];

    for (std::size_t j = 3; j < flist.size(); ++j)
      faces.push_back(flist[j]);
  }

  this->edge_info.clear();
}

/* ---------------------------------------------------------------- */

void
MeshSlicing::slice (Plane3f const& plane)
{
  this->new_verts.clear();

  MeshVertexList& verts = this->mesh->get_vertices();
  MeshFaceList& faces = this->mesh->get_faces();

  std::vector<float> dist;
  dist.resize(verts.size());

  /* Calculate distance for each point to the plane. */
  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    dist[i] = plane.point_dist(verts[i]);
    //std::cout << "Distances: " << verts[i] << " " << dist[i] << std::endl;
  }

  /* Walk over all triangles and check check for intersections. */
  std::size_t faces_size = faces.size() / 3;
  for (std::size_t i = 0; i < faces_size; ++i)
  {
    std::size_t fidx = i * 3;
    std::vector<std::pair<std::size_t, std::size_t> > intersections;
    for (std::size_t j = 0; j < 3; ++j)
    {
      std::size_t jp1 = (j + 1) % 3;
      std::size_t v1idx = faces[fidx + j];
      std::size_t v2idx = faces[fidx + jp1];
      float dist1 = dist[v1idx];
      float dist2 = dist[v2idx];

      /* Continue if both points are on the same side of the plane. */
      if ((dist1 < 0.0f && dist2 < 0.0f) || (dist1 > 0.0f && dist2 > 0.0f))
        continue;


      /* Continue if one of the points snaps to the plane. */
      #define SNAP_DIST 0.005f
      if (EPSILON_EQ(dist1, 0.0f, SNAP_DIST)
          || EPSILON_EQ(dist2, 0.0f, SNAP_DIST))
        continue;

      /* Check angle between edge and plane normal. */
      #define SNAP_ANGLE 0.2f /* acos(0.1) = 84.26. 5.74 degree eps. */
      float cosangle = plane.n * (verts[v2idx] - verts[v1idx]).norm();
      if (EPSILON_EQ(cosangle, 0.0f, SNAP_ANGLE))
        continue;

      /* Get point of intersection. */
      std::size_t vidx = this->get_intersection(v1idx, v2idx, plane);
      intersections.push_back(std::make_pair(j, vidx));
    }

    //if (!intersections.empty())
    //  std::cout << "Slicing a triangle" << std::endl;

    if (intersections.size() == 1)
    {
      std::size_t iidx = intersections.back().first;
      std::size_t v0idx = faces[fidx + iidx];
      std::size_t v1idx = intersections.back().second;
      std::size_t v2idx = faces[fidx + (iidx + 1) % 3];
      std::size_t v3idx = faces[fidx + (iidx + 2) % 3];
      faces[fidx + 0] = (MeshVIndex)v0idx;
      faces[fidx + 1] = (MeshVIndex)v1idx;
      faces[fidx + 2] = (MeshVIndex)v3idx;
      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v2idx);
      faces.push_back((MeshVIndex)v3idx);
    }
    else if (intersections.size() == 2)
    {
      std::size_t i1idx = intersections.front().first;
      std::size_t i2idx = intersections.back().first;
      std::size_t v0idx = faces[fidx + i1idx];
      std::size_t v1idx = intersections.front().second;
      std::size_t v2idx = faces[fidx + (i1idx + 1) % 3];
      std::size_t v3idx = intersections.back().second;
      std::size_t v4idx = faces[fidx + (i1idx + 2) % 3];
      if (i2idx != i1idx + 1)
        std::swap(v3idx, v4idx);

      faces[fidx + 0] = (MeshVIndex)v0idx;
      faces[fidx + 1] = (MeshVIndex)v1idx;
      faces[fidx + 2] = (MeshVIndex)v4idx;
      /* Fixme: Use best angle split */
      faces.push_back((MeshVIndex)v1idx);
      faces.push_back((MeshVIndex)v2idx);
      faces.push_back((MeshVIndex)v3idx);
      faces.push_back((MeshVIndex)v3idx);
      faces.push_back((MeshVIndex)v4idx);
      faces.push_back((MeshVIndex)v1idx);
    }
    else if (intersections.size() > 2)
      throw Exception("Too much intersections");
  }

  this->new_verts.clear();
}

/* ---------------------------------------------------------------- */

std::size_t
MeshSlicing::get_intersection (std::size_t v1idx,
    std::size_t v2idx, Plane3f const& plane)
{
  MeshVertexList& verts = this->mesh->get_vertices();

  /* Try to look up a vertex for this edge. */
  if (v1idx > v2idx)
    std::swap(v1idx, v2idx);

  std::size_t vidx;
  NewV::iterator i = this->new_verts.find(std::make_pair(v1idx, v2idx));
  if (i == this->new_verts.end())
  {
    Vec3f const& v1 = verts[v1idx];
    Vec3f d = verts[v2idx] - v1;
    float lambda = (-v1 * plane.n + plane.d) / (d * plane.n);
    Vec3f intersection = v1 + d * lambda;
    verts.push_back(intersection);
    vidx = verts.size() - 1;
    this->new_verts.insert(std::make_pair(std::make_pair(v1idx, v2idx), vidx));
  }
  else
  {
    vidx = i->second;
  }

  return vidx;
}

REMESHER_NAMESPACE_END
