#ifndef REMESHER_PLANE_HEADER
#define REMESHER_PLANE_HEADER

#include "vec3.h"
#include "defines.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Plane;

typedef Plane<float> Plane3f;
typedef Plane<double> Plane3d;
typedef Plane<int> Plane3i;
typedef Plane<unsigned int> Plane3ui;
typedef Plane<short> Plane3s;
typedef Plane<unsigned short> Plane3us;
typedef Plane<char> Plane3c;
typedef Plane<unsigned char> Plane3uc;

/* Class that represents a plane. The normal is expected to be normalized. */
template <class T>
class Plane
{
  public:
    Vec3<T> n;
    T d;

  public:
    Plane (void) {}
    /* Creates a plane with normal n and distance d from the origin. */
    Plane (Vec3<T> const& n, T const& d) : n(n), d(d) {}
    /* Creates a plane containing p with normal n. */
    Plane (Vec3<T> const& n, Vec3<T> const& p) : n(n) { d = p * n; }

    /* Returns the signed distance from a point to the plane. */
    T point_dist (Vec3<T> const& p) const
    { return p * n - d; }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_PLANE_HEADER */
