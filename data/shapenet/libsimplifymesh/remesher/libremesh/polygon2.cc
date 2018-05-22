#include "exception.h"
#include "straightline2.h"
#include "polygon2.h"

REMESHER_NAMESPACE_BEGIN

void
Polygon2::clip_with_line (Vec2f const& s1, Vec2f const& s2)
{
  /* We create a new resulting polygon during clipping. */
  Polygon2 p;

  /* Create a parameter line representation for the clipping line.
   * s1 is the base point, u is the direction vector. */
  Vec2f u = s2 - s1;
  StraightLine2f clip_line(s1, s2);

  /* Iterate over all polygon edges and create the new polygon. */
  for (std::size_t i = 0; i < this->size(); ++i)
  {
    std::size_t ip1 = (i + 1) % this->size();

    /* Create a parameter line representation for the polygon edge. */
    Vec2f const& e1 = this->at(i);
    Vec2f const& e2 = this->at(ip1);

    /* Solve edge equation to determine relative position. */
    float eeq1 = clip_line.edge_equation(e1);
    float eeq2 = clip_line.edge_equation(e2);

    /* If both polygon vertices are on the right side, skip that edge. */
    if (eeq1 > 0.0f && eeq2 > 0.0f)
      continue;

    /* If the first vertex is on the left side, append it to the polygon. */
    if (eeq1 <= 0.0f)
      p.push_back(e1);

    /* If there is a intersection with the clip line, add the intersection. */
    if ((eeq1 < 0.0f && eeq2 > 0.0f) || (eeq1 > 0.0f && eeq2 < 0.0f))
    {
      /* This is the direction of the polygon edge to be clipped. */
      Vec2f v = e2 - e1;

      /* Calculate the inverse determinante (of the 2x2 LSE). */
      float idet = u[1] * v[0] - u[0] * v[1];
      if (idet == 0.0f)
        throw Exception("Intersection of parallel lines");

      /* Calculate the intersection parameters. */
      Vec2f es = e1 - s1;
      //float lambda_s = (v[0] * es[1] - v[1] * es[0]) / idet;
      float lambda_e = (u[0] * es[1] - u[1] * es[0]) / idet;

      Vec2f intersection = e1 + v * lambda_e;
      p.push_back(intersection);
    }
  }

  (*this) = p;
}

REMESHER_NAMESPACE_END
