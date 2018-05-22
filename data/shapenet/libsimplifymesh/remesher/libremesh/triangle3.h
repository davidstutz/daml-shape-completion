#ifndef REMESHER_TRIANGLE3_HEADER
#define REMESHER_TRIANGLE3_HEADER

#include "defines.h"
#include "vec2.h"
#include "vec3.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Triangle3;
typedef Triangle3<float> Triangle3f;
typedef Triangle3<double> Triangle3d;

/* A collection of some formulars for triangles.
 * Most of these methods are quite expensive. Use them carefully!
 */
template <class T>
class Triangle3
{
  public:
    Vec3<T> a, b, c;

  public:
    Triangle3 (void);
    Triangle3 (Vec3<T> const& a, Vec3<T> const& b, Vec3<T> const& c);
    Triangle3 (TriangleMeshPtr mesh, std::size_t face);

    void set (Vec3<T> const& a, Vec3<T> const& b, Vec3<T> const& c);

    /* Returns the point for the given barycentric coordinates. */
    Vec3<T> get_point_for_bary (Vec3<T> const& bary);
    Vec3<T> get_point_for_bary (Vec2<T> const& bary);

    /* Circumcircle methods. */
    void get_circumcircle (T& qradius, Vec3<T>& center);
    Vec3<T> get_circumcenter (void);
    T get_circumradius (void);

    /* Incircle methods. */
    T get_inradius (void);

    /* Misc methods. */
    T get_quality_ratio (void);
    T get_area (void);
};

/* ---------------------------------------------------------------- */

template <class T>
Triangle3<T>::Triangle3 (void)
{
}

template <class T>
Triangle3<T>::Triangle3 (Vec3<T> const& a, Vec3<T> const& b, Vec3<T> const& c)
  : a(a), b(b), c(c)
{
}

template <class T>
Triangle3<T>::Triangle3 (TriangleMeshPtr mesh, std::size_t face)
{
  MeshVertexList const& verts = mesh->get_vertices();
  MeshFaceList const& faces = mesh->get_faces();
  Vec3f const& x = verts[faces[face * 3 + 0]];
  Vec3f const& y = verts[faces[face * 3 + 1]];
  Vec3f const& z = verts[faces[face * 3 + 2]];
  this->a[0] = (T)x[0]; this->a[1] = (T)x[1]; this->a[2] = (T)x[2];
  this->b[0] = (T)y[0]; this->b[1] = (T)y[1]; this->b[2] = (T)y[2];
  this->c[0] = (T)z[0]; this->c[1] = (T)z[1]; this->c[2] = (T)z[2];
}

template <class T>
inline void
Triangle3<T>::set (Vec3<T> const& a, Vec3<T> const& b, Vec3<T> const& c)
{
  this->a = a;
  this->b = b;
  this->c = c;
}

template <class T>
inline Vec3<T>
Triangle3<T>::get_point_for_bary (Vec3<T> const& bary)
{
  return a * bary[0] + b * bary[1] + c * bary[2];
}

template <class T>
inline Vec3<T>
Triangle3<T>::get_point_for_bary (Vec2<T> const& bary)
{
  return a * bary[0] + b * bary[1] + c * ((T)1.0 - bary[0] - bary[1]);
}

template <class T>
inline T
Triangle3<T>::get_quality_ratio (void)
{
  float l1 = (this->b - this->a).length();
  float l2 = (this->c - this->b).length();
  float l3 = (this->a - this->c).length();
  float peri2 = (l1 + l2 + l3) / 2.0f;
  float area = ::sqrtf(peri2 * (peri2 - l1) * (peri2 - l2) * (peri2 - l3));

  float circumradius = l1 * l2 * l3 / ::sqrtf
      ((l1 + l2 + l3) * (l2 + l3 - l1) * (l3 + l1 - l2) * (l1 + l2 - l3));
  float inradius = area / peri2;
  return circumradius / inradius;
}

template <class T>
inline void
Triangle3<T>::get_circumcircle (T& qradius, Vec3<T>& center)
{
  Vec3<T> u = b - a;
  Vec3<T> v = c - a;
  Vec3<T> w = u.cross(v);
  Vec3<T> c = (v * (u * u) - u * (v * v)).cross(w) / ((T)2 * (w * w));
  center = a + c;
  qradius = c * c;
}

template <class T>
inline Vec3<T>
Triangle3<T>::get_circumcenter (void)
{
  T qr;
  Vec3<T> c;
  this->get_circumcircle(qr, c);
  return c;
}

template <class T>
inline T
Triangle3<T>::get_circumradius (void)
{
  T qr;
  Vec3<T> c;
  this->get_circumcircle(qr, c);
  return std::sqrt(qr);
}

template <class T>
inline T
Triangle3<T>::get_inradius (void)
{
  float l1 = (this->b - this->a).length();
  float l2 = (this->c - this->b).length();
  float l3 = (this->a - this->c).length();
  float peri2 = (l1 + l2 + l3) / 2.0f;
  float qtarea = peri2 * (peri2 - l1) * (peri2 - l2) * (peri2 - l3);
  qtarea = MY_MAX(0.0f, qtarea);
  float area = std::sqrt(qtarea);
  return area / peri2;
}

template <class T>
inline T
Triangle3<T>::get_area (void)
{
  float l1 = (b - a).length();
  float l2 = (c - b).length();
  float l3 = (a - c).length();
  float peri2 = (l1 + l2 + l3) / 2.0f;
  float qtarea = peri2 * (peri2 - l1) * (peri2 - l2) * (peri2 - l3);
  qtarea = MY_MAX(0.0f, qtarea);
  float area = std::sqrt(qtarea);
  return area;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_TRIANGLE3_HEADER */
