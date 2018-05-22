#include <iostream>

#include "exception.h"
#include "micropatch.h"

REMESHER_NAMESPACE_BEGIN

void
MicroPatch::calculate_patch (void)
{
  /* At least four vertices are required, the center and three adjacent. */
  if (this->adj3d.size() < 3)
    throw Exception("Invalid patch with too few vertices");

  /* A list of edge length values is cached. */
  std::vector<float> edge_length;
  for (std::size_t i = 0; i < this->adj3d.size(); ++i)
    edge_length.push_back((this->adj3d[i] - this->center).length());

  /* Sum up all angles between the edges. */
  float total_angle = 0.0f;
  for (std::size_t i = 0; i < this->adj3d.size(); ++i)
  {
    std::size_t ip1 = (i + 1) % this->adj3d.size();
    Vec3f e1(this->adj3d[i] - this->center);
    Vec3f e2(this->adj3d[ip1] - this->center);
    float cosa = e1 * e2 / (edge_length[i] * edge_length[ip1]);
    //std::cout << "Angle is: " << MY_RAD2DEG(std::acos(cosa)) << std::endl;
    total_angle += std::acos(cosa);
  }

  //std::cout << "Total angle is: " << MY_RAD2DEG(total_angle) << std::endl;

  /* Create new patch vertices. */
  this->adj2d.clear();
  float angle_factor = (float)MY_2PI / total_angle;
  float current_angle = 0.0f;
  for (std::size_t i = 0; i < this->adj3d.size(); ++i)
  {
    /* Place first vertex. We create a vector (x, 0)
     * and rotate it with the current angle. */
    this->adj2d.push_back(Vec2f(edge_length[i] * std::cos(current_angle),
        edge_length[i] * std::sin(current_angle)));

    /* Calculate the angle between edges and add it to current angle. */
    std::size_t ip1 = (i + 1) % this->adj3d.size();
    Vec3f e1(this->adj3d[i] - this->center);
    Vec3f e2(this->adj3d[ip1] - this->center);
    float cosa = e1 * e2 / (edge_length[i] * edge_length[ip1]);
    float angle2d = angle_factor * std::acos(cosa);
    current_angle += angle2d;
  }
}

/* ---------------------------------------------------------------- */

TriangleMeshPtr
MicroPatch::get_debug_mesh (void) const
{
  TriangleMeshPtr mesh(TriangleMesh::create());
  MeshVertexList& verts = mesh->get_vertices();
  MeshFaceList& faces = mesh->get_faces();

#define BOTH_MICRO_PATCHES 0
#if BOTH_MICRO_PATCHES
  float xmin = 999999.0f;
  float xmax = -999999.0f;
  for (std::size_t i = 0; i < this->adj2d.size(); ++i)
  {
    xmin = (this->adj2d[i][0] < xmin ? this->adj2d[i][0] : xmin);
    xmax = (this->adj2d[i][0] > xmax ? this->adj2d[i][0] : xmax);
  }
  float xoff = (xmax - xmin) / 2.0f;

  verts.push_back(Vec3f(-xoff, 0.0f, 0.0f));
  for (std::size_t i = 0; i < this->adj3d.size(); ++i)
  {
    verts.push_back(this->adj3d[i] - this->center);
    verts.back()[0] -= xoff;
  }

  verts.push_back(Vec3f(xoff, 0.0f, 0.0f));
  for (std::size_t i = 0;  i < this->adj2d.size(); ++i)
  {
    verts.push_back(Vec3f(this->adj2d[i][0], this->adj2d[i][1], 0.0f));
    verts.back()[0] += xoff;
  }

  for (std::size_t i = 1; i < this->adj3d.size() + 1; ++i)
  {
    faces.push_back(0);
    faces.push_back(i);
    faces.push_back(1 + i % this->adj3d.size());
  }

  std::size_t ioff = this->adj3d.size() + 1;
  for (std::size_t i = 1; i < this->adj2d.size() + 1; ++i)
  {
    faces.push_back(0 + ioff);
    faces.push_back(i + ioff);
    faces.push_back(1 + i % this->adj2d.size() + ioff);
  }
#else
  verts.push_back(Vec3f(0.0f, 0.0f, 0.0f));
  for (std::size_t i = 0;  i < this->adj2d.size(); ++i)
    verts.push_back(Vec3f(this->adj2d[i][0], this->adj2d[i][1], 0.0f));

  for (std::size_t i = 1; i < this->adj2d.size() + 1; ++i)
  {
    faces.push_back(0);
    faces.push_back((MeshVIndex)i);
    faces.push_back(1 + (MeshVIndex)(i % this->adj2d.size()));
  }
#endif

  return mesh;
}

REMESHER_NAMESPACE_END
