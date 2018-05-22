#include <iostream>

#include "meshstats.h"

REMESHER_NAMESPACE_BEGIN

void
MeshStats::print_stats (void)
{
  MeshVertexList& verts = this->mesh->get_vertices();
  MeshFaceList& faces = this->mesh->get_faces();

  /* Collect information regarding vertex regularity. */
  std::size_t total_verts = 0;
  std::size_t low_degree = 0;
  std::size_t high_degree = 0;

  for (std::size_t i = 0; i < verts.size(); ++i)
  {
    switch (vinfo[i].vclass)
    {
      case VERTEX_CLASS_SIMPLE:
        if (vinfo[i].adj_faces.size() > 6)
          high_degree += 1;
        else if (vinfo[i].adj_faces.size() < 6)
          low_degree += 1;
        total_verts += 1;
        break;

      case VERTEX_CLASS_BORDER:
        if (vinfo[i].adj_faces.size() > 4)
          high_degree += 1;
        else if (vinfo[i].adj_faces.size() < 4)
          low_degree += 1;
        total_verts += 1;
        break;

      default:
        total_verts += 1;
        break;
    }
  }

  std::size_t irregular = low_degree + high_degree;
  std::size_t regular = total_verts - irregular;

  /* Collect information regarding triangle angles. */
  float min_angle = (float)MY_2PI;
  float max_angle = 0.0f;
  float avg_angle = 0.0f;
  std::size_t min_angle_tri_id = MAX_SIZE_T;
  std::size_t max_angle_tri_id = MAX_SIZE_T;

  std::size_t face_amount = faces.size() / 3;
  for (std::size_t i = 0; i < face_amount; ++i)
  {
    Vec3f v[3];
    v[0] = verts[faces[i * 3 + 0]];
    v[1] = verts[faces[i * 3 + 1]];
    v[2] = verts[faces[i * 3 + 2]];

    float smallest_angle = (float)MY_2PI;
    for (std::size_t j = 0; j < 3; ++j)
    {
      Vec3f e1 = (v[(j + 1) % 3] - v[j]).norm();
      Vec3f e2 = (v[(j + 2) % 3] - v[j]).norm();
      float scalar = e1 * e2;
      scalar = MY_MAX(-1.0f, scalar);
      scalar = MY_MIN(1.0f, scalar);
      float angle = std::acos(scalar);

      if (angle < min_angle)
      {
        min_angle = angle;
        min_angle_tri_id = i;
      }
      if (angle > max_angle)
      {
        max_angle = angle;
        max_angle_tri_id = i;
      }
      if (angle < smallest_angle)
        smallest_angle = angle;
    }

    avg_angle += smallest_angle;
  }

  avg_angle /= (float)face_amount;

  /* Print information. */
  std::cout << "Mesh statistics (" << verts.size() << " vertices, "
      << face_amount << " faces)" << std::endl
      << "  Vertices with degree: " << total_verts << std::endl
      << "  Vertices with high degree: " << high_degree << std::endl
      << "  Vertices with low degree: " << low_degree << std::endl
      << "  Irregular vertices: " << irregular
      << " (" << (100.0f * (float)irregular / (float)total_verts)
      << "%), regular vertices: " << regular
      << " (" << (100.0f * (float)regular / (float)total_verts)
      << "%)" << std::endl
      << std::endl
      << "  Minimum angle: " << MY_RAD2DEG(min_angle)
      << " (Triangle " << min_angle_tri_id << ")" << std::endl
      << "  Maximum angle: " << MY_RAD2DEG(max_angle)
      << " (Triangle " << max_angle_tri_id << ")" << std::endl
      << "  Average minimum angle: " << MY_RAD2DEG(avg_angle)
      << std::endl;
}

REMESHER_NAMESPACE_END
