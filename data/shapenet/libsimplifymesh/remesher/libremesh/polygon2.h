#ifndef REMESHER_POLYGON_2_HEADER
#define REMESHER_POLYGON_2_HEADER

#include <vector>

#include "defines.h"
#include "vec2.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class represents a polygon in two dimenstions.
 */
class Polygon2 : public std::vector<Vec2f>
{
  public:
    /* Clips the polygon along an edge specified with two points.
     * All polygon vertices on the right side of the line are removed
     * (so the ordering of clipping vertices matters), new polygon
     * vertices are introduced at the intersections. */
    void clip_with_line (Vec2f const& s1, Vec2f const& s2);
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_POLYGON_2_HEADER */
