#ifndef REMESHER_VEC3_HEADER
#define REMESHER_VEC3_HEADER

#include <ostream>
#include <cmath>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Vec3;

typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;
typedef Vec3<unsigned int> Vec3ui;
typedef Vec3<short> Vec3s;
typedef Vec3<unsigned short> Vec3us;
typedef Vec3<char> Vec3c;
typedef Vec3<unsigned char> Vec3uc;

template <class T>
class Vec3
{
  public:
    T x, y, z;

  public:
    Vec3 (void) {}
    Vec3 (T* v) : x(v[0]), y(v[1]), z(v[2]) {}
    Vec3 (T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}

    Vec3<T> plus (Vec3<T> const& c) const
    { return Vec3<T>(x + c.x, y + c.y, z + c.z); }

    Vec3<T> minus (Vec3<T> const& c) const
    { return Vec3<T>(x - c.x, y - c.y, z - c.z); }

    Vec3<T> mult (T const& f) const
    { return Vec3<T>(x * f, y * f, z * f); }

    Vec3<T> div (T const& f) const
    { return Vec3<T>(x / f, y / f, z / f); }

    Vec3<T> cross (Vec3<T> const& c) const
    { return Vec3<T>(y * c.z - z * c.y, c.x * z - x * c.z, x * c.y - c.x * y); }

    T scalar (Vec3<T> const& c) const
    { return x * c.x + y * c.y + z * c.z; }

    Vec3<T> neg (void) const
    { return Vec3<T>(-x, -y, -z); }

    T length (void) const
    { return std::sqrt(x * x + y * y + z * z); }

    T qlength (void) const
    { return x * x + y * y + z * z; }

    T* array (void)
    { return &x; }

    T const* array (void) const
    { return &x; }

    Vec3<T> norm (void) const
    { T len = length(); return Vec3<T>(x / len, y / len, z / len); }

    void norm_self (void)
    { T len = length(); x /= len; y /= len; z /= len; }

    Vec3<T> operator+ (Vec3<T> const& rhs) const
    { return plus(rhs); }

    Vec3<T>& operator+= (Vec3<T> const& rhs)
    { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }

    T const& operator[] (std::size_t i) const
    { return (&x)[i]; }

    T& operator[] (std::size_t i)
    { return (&x)[i]; }

    Vec3<T> operator- (void) const
    { return neg(); }

    Vec3<T> operator- (Vec3<T> const& rhs) const
    { return minus(rhs); }

    Vec3<T>& operator-= (Vec3<T> const& rhs)
    { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }

    Vec3<T> operator* (T const& rhs) const
    { return mult(rhs); }

    T operator* (Vec3<T> const& rhs) const
    { return scalar(rhs); }

    Vec3<T>& operator*= (T const& rhs)
    { x *= rhs; y *= rhs; z *= rhs; return *this; }

    Vec3<T> operator/ (T const& rhs) const
    { return div(rhs); }

    Vec3<T>& operator/= (T const& rhs)
    { x /= rhs; y /= rhs; z /= rhs; return *this; }

    bool operator== (Vec3<T> const& rhs)
    { return (x == rhs.x && y == rhs.y && z == rhs.z); }

    bool operator!= (Vec3<T> const& rhs)
    { return (x != rhs.x || y != rhs.y || z != rhs.z); }

    // Not compatible with STL containers
    //T* operator& (void) { return &x; }

    friend std::ostream& operator<< (std::ostream& os, Vec3<T> const& v)
    { return os << "(" << v.x << "," << v.y << "," << v.z << ")"; }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_VEC3_HEADER */
