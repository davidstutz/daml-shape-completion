#ifndef REMESHER_AVERAGE_PLANE_HEADER
#define REMESHER_AVERAGE_PLANE_HEADER

#include <vector>

#include "vec3.h"
#include "plane.h"
#include "defines.h"

/* Use GSL or GMM to solve the system. */
#define AVG_PLANE_USE_GMM 1

REMESHER_NAMESPACE_BEGIN

class AveragePlane
{
  private:
    std::vector<Vec3f> points;

  public:
    /* Copying points to this class. */
    void set_points (std::vector<Vec3f> const& points);
    void add_point (Vec3f const& point);
    void clear_points (void);

    /* Calculate average plane from copied points. */
    Plane3f get_average (void);

    /* Calculate average plane from a given set of points. */
    Plane3f get_average (std::vector<Vec3f> const& points);
};

/* ---------------------------------------------------------------- */

inline void
AveragePlane::set_points (std::vector<Vec3f> const& points)
{
  this->points = points;
}

inline void
AveragePlane::clear_points (void)
{
  this->points.clear();
}

inline void
AveragePlane::add_point (Vec3f const& point)
{
  this->points.push_back(point);
}

inline Plane3f
AveragePlane::get_average (void)
{
  return this->get_average(this->points);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_AVERAGE_PLANE_HEADER */
