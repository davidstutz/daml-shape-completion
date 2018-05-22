#ifndef REMESHER_STRAIGHT_LINE_2_HEADER
#define REMESHER_STRAIGHT_LINE_2_HEADER

#include "defines.h"
#include "vec2.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class StraightLine2;
typedef StraightLine2<float> StraightLine2f;
typedef StraightLine2<double> StraightLine2d;
typedef StraightLine2<int> StraightLine2i;
typedef StraightLine2<unsigned int> StraightLine2ui;
typedef StraightLine2<short> StraightLine2s;
typedef StraightLine2<unsigned short> StraightLine2us;
typedef StraightLine2<char> StraightLine2c;
typedef StraightLine2<unsigned char> StraightLine2uc;

template <class T>
class StraightLine2
{
  private:
    Vec2<T> a, b;

  public:
    StraightLine2 (void);
    StraightLine2 (Vec2<T> const& a, Vec2<T> const& b);

    /* The edge equation. Returns a positive value if (x,y) is right of
     * line a -> b, a nevative value if (x,y) is left of the line, and
     * zero if (x,y) is on the line. */
    T edge_equation (T const& x, T const& y);
    T edge_equation (Vec2<T> const& p);

    /* Calculates the square point distance to the line. */
    T point_qdist (Vec2<T> const& p);
};

/* ---------------------------------------------------------------- */

template <class T>
inline
StraightLine2<T>::StraightLine2 (void)
{
}

template <class T>
inline
StraightLine2<T>::StraightLine2 (Vec2<T> const& a, Vec2<T> const& b)
  : a(a), b(b)
{
}

template <class T>
inline T
StraightLine2<T>::edge_equation (T const& x, T const& y)
{
  return (T)((x - a[0]) * (b[1] - a[1]) - (y - a[1]) * (b[0] - a[0]));
}

template <class T>
inline T
StraightLine2<T>::edge_equation (Vec2<T> const& p)
{
  return this->edge_equation(p[0], p[1]);
}

template <class T>
inline T
StraightLine2<T>::point_qdist (Vec2<T> const& p)
{
  Vec2<T> r = b - a;
  T lambda = r.scalar(p - a) / r.scalar(r);
  Vec2<T> c = a + r * lambda;
  return (p - c).qlength();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_STRAIGHT_LINE_2_HEADER */
