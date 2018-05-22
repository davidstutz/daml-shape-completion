#ifndef REMESHER_TRIANGLE2_HEADER
#define REMESHER_TRIANGLE2_HEADER

#include "defines.h"
#include "vec2.h"
#include "vec3.h"
#include "matrix2.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Triangle2;
typedef Triangle2<float> Triangle2f;
typedef Triangle2<double> Triangle2d;

/* A collection of some formulars for triangles.
 * Most of these methods are quite expensive. Use them carefully!
 */
template <class T>
class Triangle2
{
  public:
    Vec2<T> a, b, c;

  public:
    Triangle2 (void);
    Triangle2 (Vec2<T> const& a, Vec2<T> const& b, Vec2<T> const& c);
    Triangle2 (TriangleMeshPtr mesh, std::size_t face);

    void set (Vec2<T> const& a, Vec2<T> const& b, Vec2<T> const& c);

    /* Returns the area of the triangle. */
    T get_area (void);

    /* Returns the triangles barycentric coordinates of x. */
    Vec3<T> get_bary_coords (Vec2<T> const& x) const;

    /* Returns the point for the given barycentric coordinates. */
    Vec2<T> get_point_for_bary (Vec3<T> const& b) const;
    Vec2<T> get_point_for_bary (Vec2<T> const& b) const;

    /* Returns the quadratic radius and center of the circumradius.
     * This function is efficient and does not use any roots. */
    void get_circumcircle (T& qradius, Vec2<T>& center) const;

    /* Implementation using the get_circumcircle() method above.
     * Hopefully, the copiler removes the small overhead. */
    Vec2<T> get_circumcenter (void) const;
    T get_circumradius (void) const;
};

/* ---------------------------------------------------------------- */

template <class T>
inline
Triangle2<T>::Triangle2 (void)
{
}

template <class T>
inline
Triangle2<T>::Triangle2 (Vec2<T> const& a, Vec2<T> const& b, Vec2<T> const& c)
  : a(a), b(b), c(c)
{
}

template <class T>
Triangle2<T>::Triangle2 (TriangleMeshPtr mesh, std::size_t face)
{
  MeshVertexList const& verts = mesh->get_vertices();
  MeshFaceList const& faces = mesh->get_faces();
  Vec3f const& x = verts[faces[face * 3 + 0]];
  Vec3f const& y = verts[faces[face * 3 + 1]];
  Vec3f const& z = verts[faces[face * 3 + 2]];
  this->a[0] = (T)x[0]; this->a[1] = (T)x[1];
  this->b[0] = (T)y[0]; this->b[1] = (T)y[1];
  this->c[0] = (T)z[0]; this->c[1] = (T)z[1];
}

template <class T>
inline T
Triangle2<T>::get_area (void)
{
  return (b[0] * c[1] + c[0] * a[1] + a[0] * b[1]
      - a[1] * b[0] - b[1] * c[0] - c[1] * a[0])/* / (T)2*/;
}

template <class T>
inline void
Triangle2<T>::set (Vec2<T> const& a, Vec2<T> const& b, Vec2<T> const& c)
{
  this->a = a;
  this->b = b;
  this->c = c;
}

template <class T>
inline Vec3<T>
Triangle2<T>::get_bary_coords (Vec2<T> const& x) const
{
  Matrix2<T> m;
  m[0] = a[0] - c[0];  m[1] = b[0] - c[0];
  m[2] = a[1] - c[1];  m[3] = b[1] - c[1];
  m = m.invert();
  Vec2<T> bary = m * (x - c);
  return Vec3<T>(bary[0], bary[1], (T)1 - bary[0] - bary[1]);
}

template <class T>
inline Vec2<T>
Triangle2<T>::get_point_for_bary (Vec3<T> const& bc) const
{
  return a * bc[0] + b * bc[1] + c * bc[2];
}

template <class T>
inline Vec2<T>
Triangle2<T>::get_point_for_bary (Vec2<T> const& bc) const
{
  return a * bc[0] + b * bc[1] + c * ((T)1.0 - bc[0] - bc[1]);
}

template <class T>
inline void
Triangle2<T>::get_circumcircle (T& qradius, Vec2<T>& center) const
{
  T qa = a[0] * a[0] + a[1] * a[1];
  T qb = b[0] * b[0] + b[1] * b[1];
  T qc = c[0] * c[0] + c[1] * c[1];

  T sx = (qa * b[1] + qb * c[1] + qc * a[1]
      - qc * b[1] - qb * a[1] - qa * c[1]) / (T)2;
  T sy = (a[0] * qb + b[0] * qc + c[0] * qa
      - c[0] * qb - b[0] * qa - a[0] * qc) / (T)2;

  T sa = a[0] * b[1] + b[0] * c[1] + c[0] * a[1]
      - c[0] * b[1] - b[0] * a[1] - a[0] * c[1];

  center[0] = sx / sa;
  center[1] = sy / sa;

  T sb = a[0] * b[1] * qc + b[0] * c[1] * qa + c[0] * a[1] * qb
      - c[0] * b[1] * qa - b[0] * a[1] * qc - a[0] * c[1] * qb;

  qradius = sb / sa + (sx * sx + sy * sy) / (sa * sa);
}

template <class T>
inline Vec2<T>
Triangle2<T>::get_circumcenter (void) const
{
  float qradius;
  Vec2<T> center;
  this->get_circumcircle(qradius, center);
  return center;
}

template <class T>
inline T
Triangle2<T>::get_circumradius (void) const
{
  float qradius;
  Vec2<T> center;
  this->get_circumcircle(qradius, center);
  return std::sqrt(qradius);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_TRIANGLE2_HEADER */
