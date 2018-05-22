#ifndef REMESHER_VEC2_HEADER
#define REMESHER_VEC2_HEADER

#include <ostream>
#include <cmath>
#include <cfloat>

#include "defines.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Vec2;

typedef Vec2<float> Vec2f;
typedef Vec2<double> Vec2d;
typedef Vec2<int> Vec2i;
typedef Vec2<unsigned int> Vec2ui;
typedef Vec2<short> Vec2s;
typedef Vec2<unsigned short> Vec2us;
typedef Vec2<char> Vec2c;
typedef Vec2<unsigned char> Vec2uc;

template <class T>
class Vec2
{
  public:
    T x, y;

  public:
    Vec2 (void) {}
    Vec2 (T* v) : x(v[0]), y(v[1]) {}
    Vec2 (T _x, T _y) : x(_x), y(_y) {}

    Vec2<T> plus (Vec2<T> const& c) const
    { return Vec2<T>(x + c.x, y + c.y); }

    Vec2<T> minus (Vec2<T> const& c) const
    { return Vec2<T>(x - c.x, y - c.y); }

    Vec2<T> mult (T const& f) const
    { return Vec2<T>(x * f, y * f); }

    Vec2<T> div (T const& f) const
    { return Vec2<T>(x / f, y / f); }

    T scalar (Vec2<T> const& c) const
    { return x * c.x + y * c.y; }

    Vec2<T> neg (void) const
    { return Vec2<T>(-x, -y); }

    T length (void) const
    { return std::sqrt(x * x + y * y); }

    T qlength (void) const
    { return x * x + y * y; }

    T* array (void)
    { return &x; }

    T const* array (void) const
    { return &x; }

    Vec2<T> norm (void) const
    { T len = length(); return Vec2<T>(x / len, y / len); }

    void norm_self (void)
    { T len = length(); x /= len; y /= len; }

    bool eps_equals (Vec2<T> const& c) const;

    Vec2<T> operator+ (Vec2<T> const& rhs) const
    { return plus(rhs); }

    Vec2<T>& operator+= (Vec2<T> const& rhs)
    { x += rhs.x; y += rhs.y; return *this; }

    T const& operator[] (std::size_t i) const
    { return (&x)[i]; }

    T& operator[] (std::size_t i)
    { return (&x)[i]; }

    Vec2<T> operator- (void) const
    { return neg(); }

    Vec2<T> operator- (Vec2<T> const& rhs) const
    { return minus(rhs); }

    Vec2<T>& operator-= (Vec2<T> const& rhs)
    { x -= rhs.x; y -= rhs.y; return *this; }

    Vec2<T> operator* (T const& rhs) const
    { return mult(rhs); }

    T operator* (Vec2<T> const& rhs) const
    { return scalar(rhs); }

    Vec2<T>& operator*= (T const& rhs)
    { x *= rhs; y *= rhs; return *this; }

    Vec2<T> operator/ (T const& rhs) const
    { return div(rhs); }

    Vec2<T>& operator/= (T const& rhs)
    { x /= rhs; y /= rhs; return *this; }

    // Not compatible with STL containers
    //T* operator& (void) { return &x; }

    friend std::ostream& operator<< (std::ostream& os, Vec2<T> const& v)
    { return os << "(" << v.x << "," << v.y << ")"; }
};

/* ---------------------------------------------------------------- */

template <class T>
inline bool
Vec2<T>::eps_equals (Vec2<T> const& c) const
{
  return x == c.x && y == c.y;
}

template <>
inline bool
Vec2<float>::eps_equals (Vec2<float> const& c) const
{
  return FLOAT_EQ(x, c.x) && FLOAT_EQ(y, c.y);
}

template <>
inline bool
Vec2<double>::eps_equals (Vec2<double> const& c) const
{
  return DOUBLE_EQ(x, c.x) && DOUBLE_EQ(y, c.y);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_VEC2_HEADER */
