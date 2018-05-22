#ifndef REMESHER_QUATERNION_HEADER
#define REMESHER_QUATERNION_HEADER

#include <cmath>
#include "defines.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Quaternion;

typedef Quaternion<float> Quat4f;
typedef Quaternion<double> Quat4d;
typedef Quaternion<int> Quat4i;
typedef Quaternion<unsigned int> Quat4ui;
typedef Quaternion<short> Quat4s;
typedef Quaternion<unsigned short> Quat4us;
typedef Quaternion<char> Quat4c;
typedef Quaternion<unsigned char> Quat4uc;

template <class T>
class Quaternion
{
  public:
    T r, x, y, z;

    Quaternion (void) {}
    Quaternion (T r, T* v) : r(r), x(v[0]), y(v[1]), z(v[2]) {}
    Quaternion (T* v) : r(v[0]), x(v[1]), y(v[2]), z(v[3]) {}
    Quaternion (T r, T x, T y, T z) : r(r), x(x), y(y), z(z) {}

    void set_axis_angle (T const* axis, T const& angle)
    {
      r = (T)::cos(angle / (T)2);
      T sa = (T)::sin(angle / (T)2);
      x = axis[0] * sa;
      y = axis[1] * sa;
      z = axis[2] * sa;
    }

    void get_axis_angle (T* axis, T& angle)
    {
      T scale = (T)::sqrt(x * x + y * y + z * z);
      if (scale == (T)0)
      {
        axis[0] = (T)1;
        axis[1] = (T)0;
        axis[2] = (T)0;
        angle = (T)0;
      }
      else
      {
        axis[0] = x / scale;
        axis[1] = y / scale;
        axis[2] = z / scale;
        angle = (T)(2.0 * ::acos(r));
      }
    }

    void rotate_vector (T const* vec, T* ret)
    {
      float xxzz = x * x - z * z;
      float rryy = r * r - y * y;
      float yyrrxxzz = y * y + r * r - x * x - z * z;

      float xr2 = x * r * (T)2;
      float xy2 = x * y * (T)2;
      float xz2 = x * z * (T)2;
      float yr2 = y * r * (T)2;
      float yz2 = y * z * (T)2;
      float zr2 = z * r * (T)2;

      ret[0] = (xxzz+rryy) * vec[0] + (xy2+zr2) * vec[1] + (xz2-yr2) * vec[2];
      ret[1] = (xy2-zr2) * vec[0] + (yyrrxxzz) * vec[1] + (yz2+xr2) * vec[2];
      ret[2] = (xz2+yr2) * vec[0] + (yz2-xr2) * vec[1] + (rryy-xxzz) * vec[2];
    }

    Quaternion<T> plus (Quaternion<T> const& c) const
    { return Quaternion<T>(r + c.r, x + c.x, y + c.y, z + c.z); }

    Quaternion<T> plus (T const& f) const
    { return Quaternion<T>(r + f, x, y, z); }

    Quaternion<T> minus (Quaternion<T> const& c) const
    { return Quaternion<T>(r - c.r, x - c.x, y - c.y, z - c.z); }

    Quaternion<T> minus (T const& f) const
    { return Quaternion<T>(r - f, x, y, z); }

    Quaternion<T> mult (Quaternion<T> const& c) const
    {
      return Quaternion<T>(
          r * c.r - x * c.x - y * c.y - z * c.z,
          r * c.x + x * c.r + y * c.z - z * c.y,
          r * c.y - x * c.z + y * c.r + z * c.x,
          r * c.z + x * c.y - y * c.x + z * c.r);
    }

    Quaternion<T> mult (T const& f)
    { return Quaternion<T>(r * f, x * f, y * f, z * f); }

    Quaternion<T> div (T const& f)
    { return Quaternion<T>(r / f, x / f, y / f, z / f); }

    Quaternion<T> conjugate (void) const
    { return Quaternion<T>(r, -x, -y, -z); }

    /* Length operations. */
    T length (void) const
    { return (T)sqrt(r * r + x * x + y * y + z * z); }

    T qlength (void) const
    { return r * r + x * x + y * y + z * z; }

    Quaternion<T> norm (void) const
    { T l = length(); return Quaternion<T>(r / l, x / l, y / l, z / l); }

    /* Array access for members. */
    T* array (void)
    { return &r; }

    T const* array (void) const
    { return &r; }

    /* Access operators. */
    T const& operator[] (std::size_t i) const
    { return (&r)[i]; }

    T& operator[] (std::size_t i)
    { return (&r)[i]; }

    /* Math operations operators. */
    Quaternion<T> operator+ (Quaternion<T> const& rhs) const
    { return plus(rhs); }

    Quaternion<T> operator+ (T const& rhs) const
    { return plus(rhs); }

    Quaternion<T> operator- (void) const
    { return conjugate(); }

    Quaternion<T> operator- (Quaternion<T> const& rhs) const
    { return minus(rhs); }

    Quaternion<T> operator- (T const& rhs) const
    { return minus(rhs); }

    Quaternion<T> operator* (Quaternion<T> const& rhs) const
    { return mult(rhs); }

    Quaternion<T> operator* (T const& rhs) const
    { return mult(rhs); }

    Quaternion<T> operator/ (T const& rhs) const
    { return div(rhs); }

    /* Operators on self. */
    Quaternion<T>& operator+= (Quaternion<T> const& rhs)
    { *this = plus(rhs); return *this; }

    Quaternion<T>& operator+= (T const& rhs)
    { *this = plus(rhs); return *this; }

    Quaternion<T>& operator-= (Quaternion<T> const& rhs)
    { *this = minus(rhs); return *this; }

    Quaternion<T>& operator-= (T const& rhs)
    { *this = minus(rhs); return *this; }

    Quaternion<T>& operator*= (Quaternion<T> const& rhs)
    { *this = mult(rhs); return *this; }

    Quaternion<T>& operator*= (T const& rhs)
    { *this = mult(rhs); return *this; }

    Quaternion<T> operator/= (T const& rhs)
    { r /= rhs.r; x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }

    /* Comparison operators. */
    Quaternion<T> operator== (Quaternion<T> const& rhs)
    { return r == rhs.r && x == rhs.x && y == rhs.y && z == rhs.z; }

    Quaternion<T> operator!= (Quaternion<T> const& rhs)
    { return !(*this == rhs); }

    friend std::ostream& operator<< (std::ostream& s, Quaternion<T> const& q)
    { return s << "(" << q.r << "|" << q.x << "," << q.y << "," << q.z << ")"; }

};

REMESHER_NAMESPACE_END

#endif /* REMESHER_QUATERNION_HEADER */

