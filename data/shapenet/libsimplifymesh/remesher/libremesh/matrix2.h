#ifndef REMESHER_MATRIX2_HEADER
#define REMESHER_MATRIX2_HEADER

#include <ostream>

#include "defines.h"
#include "vec2.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Matrix2;

typedef Matrix2<float> Matrix2f;
typedef Matrix2<double> Matrix2d;
typedef Matrix2<int> Matrix2i;
typedef Matrix2<unsigned int> Matrix2ui;
typedef Matrix2<char> Matrix2c;

template <class T>
class Matrix2
{
  public:
    T m[4];

    Matrix2 (void) {}
    Matrix2 (T* matrix)
    { for (std::size_t i = 0; i < 4; ++i) m[i] = matrix[i]; }

    Matrix2<T> plus (Matrix2<T> const& rhs) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i) ret[i] = m[i] + rhs[i];
      return ret;
    }

    Matrix2<T> minus (Matrix2<T> const& rhs) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i) ret[i] = m[i] - rhs[i];
      return ret;
    }

    Matrix2<T> mult (Matrix2<T> const& rhs) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i)
      {
        ret[i] = (T)0;
        for (std::size_t j = 0; j < 2; ++j)
          ret[i] += m[i-i%2 + j] * rhs[i%2 + j*2];
      }
      return ret;
    }

    Matrix2<T> mult (T const& rhs) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i) ret[i] = m[i] * rhs;
      return ret;
    }

    Vec2<T> mult (Vec2<T> const& rhs) const
    {
      return Vec2<T>(m[0] * rhs.x + m[1] * rhs.y, m[2] * rhs.x + m[3] * rhs.y);
    }

    Matrix2<T> div (T const& rhs) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i) ret[i] = m[i] / rhs;
      return ret;
    }

    Matrix2<T> neg (void) const
    {
      Matrix2<T> ret;
      for (std::size_t i = 0; i < 4; ++i) ret[i] = -m[i];
      return ret;
    }

    Matrix2<T> invert (void) const
    {
      Matrix2<T> ret;
      T det = m[0] * m[3] - m[1] * m[2];
      ret[0] = m[3] / det;
      ret[1] = -m[1] / det;
      ret[2] = -m[2] / det;
      ret[3] = m[0] / det;
      return ret;
    }

    T const* array (void) const
    { return m; }

    T* array (void)
    { return m; }

    Matrix2<T> operator+ (Matrix2<T> const& rhs) const
    { return plus(rhs); }

    Matrix2<T> operator- (Matrix2<T> const& rhs) const
    { return minus(rhs); }

    Matrix2<T> operator* (Matrix2<T> const& rhs) const
    { return mult(rhs); }

    Matrix2<T> operator* (T const& rhs) const
    { return mult(rhs); }

    Vec2<T> operator* (Vec2<T> const& rhs) const
    { return mult(rhs); }

    T const& operator[] (std::size_t i) const
    { return m[i]; }

    T& operator[] (std::size_t i)
    { return m[i]; }

    Matrix2<T>& operator= (Matrix2<T> const& rhs)
    {
      for (std::size_t i = 0; i < 4; ++i) m[i] = rhs[i];
      return *this;
    }

    friend std::ostream& operator<< (std::ostream& os, Matrix2<T> const& x)
    {
      return os << "[" << x[0] << ", " << x[1]
          << "; " << x[2] << ", " << x[3]
          << "]";
    }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_MATRIX2_HEADER */
