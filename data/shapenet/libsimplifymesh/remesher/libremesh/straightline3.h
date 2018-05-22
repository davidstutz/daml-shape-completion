#ifndef REMESHER_STRAIGHT_LINE_3_HEADER
#define REMESHER_STRAIGHT_LINE_3_HEADER

#include "defines.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class StraightLine3;
typedef StraightLine3<float> StraightLine3f;
typedef StraightLine3<double> StraightLine3d;
typedef StraightLine3<int> StraightLine3i;
typedef StraightLine3<unsigned int> StraightLine3ui;
typedef StraightLine3<short> StraightLine3s;
typedef StraightLine3<unsigned short> StraightLine3us;
typedef StraightLine3<char> StraightLine3c;
typedef StraightLine3<unsigned char> StraightLine3uc;

template <class T>
class StraightLine3
{
  private:
    Vec3<T> a, b;

  public:
    StraightLine3 (void);
    StraightLine3 (Vec3<T> const& a, Vec3<T> const& b);

    T point_dist (Vec3<T> const& p) const;
};

/* ---------------------------------------------------------------- */

template <class T>
inline
StraightLine3<T>::StraightLine3 (void)
{
}

template <class T>
inline
StraightLine3<T>::StraightLine3 (Vec3<T> const& a, Vec3<T> const& b)
  : a(a), b(b)
{
}

template <class T>
inline T
StraightLine3<T>::point_dist (Vec3<T> const& p) const
{
  return (b - a).cross(a - p).length() / (b - a).length();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_STRAIGHT_LINE_3_HEADER */
