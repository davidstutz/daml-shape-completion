#ifndef REMESHER_VEC4_HEADER
#define REMESHER_VEC4_HEADER

#include <ostream>
#include <cmath>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Vec4;

typedef Vec4<float> Vec4f;
typedef Vec4<double> Vec4d;
typedef Vec4<int> Vec4i;
typedef Vec4<unsigned int> Vec4ui;
typedef Vec4<short> Vec4s;
typedef Vec4<unsigned short> Vec4us;
typedef Vec4<char> Vec4c;
typedef Vec4<unsigned char> Vec4uc;

template <class T>
class Vec4
{
  public:
    T x, y, z, w;

  public:
    Vec4 (void) {}
    Vec4 (T* v) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
    Vec4 (T _x, T _y, T _z, T _w) : x(_x), y(_y), z(_z), w(_w) {}

    Vec4<T> plus (Vec4<T> const& c) const
    { return Vec4<T>(x + c.x, y + c.y, z + c.z, w + c.w); }

    Vec4<T> minus (Vec4<T> const& c) const
    { return Vec4<T>(x - c.x, y - c.y, z - c.z, w - c.w); }

    Vec4<T> mult (T const& f) const
    { return Vec4<T>(x * f, y * f, z * f, w * f); }

    Vec4<T> div (T const& f) const
    { return Vec4<T>(x / f, y / f, z / f, w / f); }

    T scalar (Vec4<T> const& c) const
    { return x * c.x + y * c.y + z * c.z + w * c.w; }

    Vec4<T> neg (void) const
    { return Vec4<T>(-x, -y, -z, -w); }

    T length (void) const
    { return std::sqrt(x * x + y * y + z * z + w * w); }

    T* array (void)
    { return &x; }

    T const* array (void) const
    { return &x; }

    Vec4<T> norm (void) const
    { T len = length(); return Vec4<T>(x / len, y / len, z / len, w / len); }

    void norm_self (void)
    { T len = length(); x /= len; y /= len; z /= len; w /= len; }

    Vec4<T> operator+ (Vec4<T> const& rhs) const
    { return plus(rhs); }

    Vec4<T>& operator+= (Vec4<T> const& rhs)
    { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }

    T const& operator[] (std::size_t i) const
    { return (&x)[i]; }

    T& operator[] (std::size_t i)
    { return (&x)[i]; }

    Vec4<T> operator- (void) const
    { return neg(); }

    Vec4<T> operator- (Vec4<T> const& rhs) const
    { return minus(rhs); }

    Vec4<T>& operator-= (Vec4<T> const& rhs)
    { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }

    Vec4<T> operator* (T const& rhs) const
    { return mult(rhs); }

    T operator* (Vec4<T> const& rhs) const
    { return scalar(rhs); }

    Vec4<T>& operator*= (T const& rhs)
    { x *= rhs; y *= rhs; z *= rhs; w *= rhs; return *this; }

    Vec4<T> operator/ (T const& rhs) const
    { return div(rhs); }

    Vec4<T>& operator/= (T const& rhs)
    { x /= rhs; y /= rhs; z /= rhs; w /= rhs; return *this; }

    // Not compatible with STL containers
    //T* operator& (void) { return &x; }

    friend std::ostream& operator<< (std::ostream& os, Vec4<T> const& v)
    { return os << "(" << v.x << "," << v.y << "," << v.z << "," << v.w << ")"; }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_VEC4_HEADER */
