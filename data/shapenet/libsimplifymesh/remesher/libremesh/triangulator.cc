#include "plane.h"
#include "triangulator.h"

REMESHER_NAMESPACE_BEGIN

bool
Triangulator::triangulate (VertexInfo::VertexList const& vlist,
    MeshFaceList& flist, Vec3f const& ap_normal, SplitLine const& sline)
{
  /* If we have less than three vertices, we can't triangulate.
   * This happens if a "faulty" splitline is specified. */
  if (vlist.size() < 3)
    return true;

  /* If we have three vertices, we're done with the given loop. */
  if (vlist.size() == 3)
  {
    for (std::size_t i = 0; i < 3; ++i)
      flist.push_back((MeshVIndex)vlist[i]);
    return true;
  }

  float best_quality = -1;
  SplitLine sl = sline;
  std::size_t vsize = vlist.size();

  /* If a split line is specified, test if it's valid. */
  if (sline.first != sline.second)
  {
    best_quality = this->check_splitline_quality(sline, vlist, ap_normal);
  }
  else
  {
    /* If no split line is specified, choose a splitline.
     * Bulid every possible split plane and measure quality.*/
    for (std::size_t i = 0; i < vsize - 2; ++i)
      for (std::size_t j = i + 2;
          j != (i + vsize - 1) % vsize;
          j = (j + 1) % vsize)
      {
        /* Create split plane from vertex i to j. */
        float quality = this->check_splitline_quality
            (SplitLine(i, j), vlist, ap_normal);

        //std::cout << "Quality from " << vlist[i]
        //    << " to " << vlist[j] << ": " << quality <<  std::endl;
        if (quality > best_quality)
        {
          best_quality = quality;
          sl = SplitLine(i, j);
        }
      }
  }

  /* Check if we have a valid split plane. */
  if (best_quality < 0.0f)
    return false;

  //std::cout << "Best quality split plane is from " << vlist[sp_a] << " to "
  //    << vlist[sp_b] << " with quality " << best_quality << std::endl;

  /* The split plane is good! Divide and conquer by creating two
   * loops and recursively triangulate. */

  VertexInfo::VertexList loop1;
  for (std::size_t i = sl.first; i != sl.second; i = (i + 1) % vsize)
    loop1.push_back(vlist[i]);
  loop1.push_back(vlist[sl.second]);

  VertexInfo::VertexList loop2;
  for (std::size_t i = sl.second; i != sl.first; i = (i + 1) % vsize)
    loop2.push_back(vlist[i]);
  loop2.push_back(vlist[sl.first]);

  if (this->triangulate(loop1, flist, ap_normal, SplitLine(0, 0))
      && this->triangulate(loop2, flist, ap_normal, SplitLine(0, 0)))
    return true;

  return false;
}

/* ---------------------------------------------------------------- */

float
Triangulator::check_splitline_quality (SplitLine const& sl,
    VertexInfo::VertexList const& vl, Vec3f const& ap_normal)
{
  MeshVertexList const& verts = this->mesh->get_vertices();

  if (this->vinfo.get() != 0
      && this->vinfo->is_mesh_edge(vl[sl.first], vl[sl.second]))
    return -1.0f;

  Vec3f const& a(verts[vl[sl.first]]);
  Vec3f const& b(verts[vl[sl.second]]);

  int direction = 0;
  float min_distance = -1.0f;
  Plane3f split_plane((b - a).cross(ap_normal).norm(), a);

  /* Walk over all points but skip split line points and
   * measure distance. */
  for (std::size_t k = (sl.first + 1) % vl.size();
      k != sl.first; k = (k + 1) % vl.size())
  {
    /* Skip direction when passing the second split plane vertex. */
    if (k == sl.second)
    {
      direction *= -1;
      continue;
    }

    /* Measure distance. */
    float point_dist = split_plane.point_dist(verts[vl[k]]);

    /* Check if point is on the right side. Break loop otherwise. */
    if (direction == 0)
      direction = (point_dist < 0.0f ? -1 : 1);
    //else if (FLOAT_EQ(point_dist, 0.0f)
    //{
    //  return -1.0f;
    //}
    else if (direction != (point_dist < 0.0f ? -1 : 1))
    {
      /* Indicate error and return. */
      return -1.0f;
    }

    if (min_distance < 0.0f || std::fabs(point_dist) < min_distance)
      min_distance = std::fabs(point_dist);
  }

  /* Indicate error if there was not a single valid split line. */
  if (min_distance < 0.0f)
    return -1.0f;

  /* The quality for the split plane. */
  return (min_distance * min_distance) / (b - a).qlength();
}

REMESHER_NAMESPACE_END
