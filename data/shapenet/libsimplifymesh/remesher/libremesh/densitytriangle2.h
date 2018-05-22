#ifndef REMESHER_DENSITY_TRIANGLE_2_HEADER
#define REMESHER_DENSITY_TRIANGLE_2_HEADER

#include "defines.h"
#include "triangle2.h"
#include "vec2.h"
#include "vec3.h"
#include "matrix3.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class DensityTriangle2;
typedef DensityTriangle2<float> DensityTriangle2f;
typedef DensityTriangle2<double> DensityTriangle2d;

template <class T>
class DensityTriangle2 : public Triangle2<T>
{
  public:
    Vec3<T> dens;

  public:
    DensityTriangle2 (void);
    DensityTriangle2 (Vec2<T> const& a, Vec2<T> const& b,
        Vec2<T> const& c, Vec3<T> const& dens);

    void set_dens (Vec3<T> const& dens);

    /* Returns the mass of the triangle, which is the non-uniform
     * area of the triangle using the density values at the vertices. */
    T get_mass (void);

    Vec2<T> get_centroid (void);
};

/* ---------------------------------------------------------------- */

template <class T>
inline
DensityTriangle2<T>::DensityTriangle2 (void)
{
}

template <class T>
inline
DensityTriangle2<T>::DensityTriangle2 (Vec2<T> const& a, Vec2<T> const& b,
        Vec2<T> const& c, Vec3<T> const& dens)
  : Triangle2<T>(a, b, c), dens(dens)
{
}

template <class T>
inline void
DensityTriangle2<T>::set_dens (Vec3<T> const& dens)
{
  this->dens = dens;
}

template <class T>
inline T
DensityTriangle2<T>::get_mass (void)
{
  T area = this->get_area();
  return area * (dens[0] + dens[1] + dens[2]) / (T)3;
}

template <class T>
inline Vec2<T>
DensityTriangle2<T>::get_centroid (void)
{
  T area = this->get_area();
  T mass3 = area * (dens[0] + dens[1] + dens[2]);

  Vec2<T> const& a = this->a;
  Vec2<T> const& b = this->b;
  Vec2<T> const& c = this->c;

  T cx = dens[0] * a[0] + dens[1] * b[0] + dens[2] * c[0]
      + (dens[0] * b[0] + dens[1] * a[0]
      + dens[0] * c[0] + dens[2] * a[0]
      + dens[1] * c[0] + dens[2] * b[0]) / 2.0f;
  T cy = dens[0] * a[1] + dens[1] * b[1] + dens[2] * c[1]
      + (dens[0] * b[1] + dens[1] * a[1]
      + dens[0] * c[1] + dens[2] * a[1]
      + dens[1] * c[1] + dens[2] * b[1]) / 2.0f;

  return Vec2<T>(cx * area / ((T)2 * mass3), cy * area / ((T)2 * mass3));
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_DENSITY_TRIANGLE_2_HEADER */
