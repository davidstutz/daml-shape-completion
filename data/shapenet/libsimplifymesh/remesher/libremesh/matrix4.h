#ifndef REMESHER_MATRIX4_HEADER
#define REMESHER_MATRIX4_HEADER

#include <ostream>

#include "defines.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Matrix4;

typedef Matrix4<float> Matrix4f;
typedef Matrix4<double> Matrix4d;
typedef Matrix4<int> Matrix4i;
typedef Matrix4<unsigned int> Matrix4ui;
typedef Matrix4<char> Matrix4c;

template <class T>
class Matrix4
{
  public:
    T m[16];

    Matrix4 (void) {}
    Matrix4 (T* matrix)
    { for (std::size_t i = 0; i < 16; ++i) m[i] = matrix[i]; }

    Matrix4<T> plus (Matrix4<T> const& rhs) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i) ret[i] = m[i] + rhs[i];
      return ret;
    }

    Matrix4<T> minus (Matrix4<T> const& rhs) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i) ret[i] = m[i] - rhs[i];
      return ret;
    }

    Matrix4<T> mult (Matrix4<T> const& rhs) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i)
      {
        ret[i] = (T)0;
        for (std::size_t j = 0; j < 4; ++j)
          ret[i] += m[i-i%4 + j] * rhs[i%4 + j*4];
      }
      return ret;
    }

    Matrix4<T> mult (T const& rhs) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i) ret[i] = m[i] * rhs;
      return ret;
    }

    Vec3<T> mult (Vec3<T> const& rhs) const
    {
      return Vec3<T>(m[0] * rhs.x + m[1] * rhs.y + m[2] * rhs.z + m[3],
          m[4] * rhs.x + m[5] * rhs.y + m[6] * rhs.z + m[7],
          m[8] * rhs.x + m[9] * rhs.y + m[10] * rhs.z + m[11]);
    }

    Matrix4<T> div (T const& rhs) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i) ret[i] = m[i] / rhs;
      return ret;
    }

    Matrix4<T> neg (void) const
    {
      Matrix4<T> ret;
      for (std::size_t i = 0; i < 16; ++i) ret[i] = -m[i];
      return ret;
    }

    Matrix4<T> transpose (void) const
    {
      Matrix4<T> ret;
      ret[0]  = m[0]; ret[1]  = m[4]; ret[2]  = m[8];  ret[3]  = m[12];
      ret[4]  = m[1]; ret[5]  = m[5]; ret[6]  = m[9];  ret[7]  = m[13];
      ret[8]  = m[2]; ret[9]  = m[6]; ret[10] = m[10]; ret[11] = m[14];
      ret[12] = m[3]; ret[13] = m[7]; ret[14] = m[11]; ret[15] = m[15];
      return ret;
    }

    T const* array (void) const
    { return m; }

    T* array (void)
    { return m; }

    Matrix4<T> operator+ (Matrix4<T> const& rhs) const
    { return plus(rhs); }

    Matrix4<T> operator- (Matrix4<T> const& rhs) const
    { return minus(rhs); }

    Matrix4<T> operator* (Matrix4<T> const& rhs) const
    { return mult(rhs); }

    Matrix4<T> operator* (T const& rhs) const
    { return mult(rhs); }

    Vec3<T> operator* (Vec3<T> const& rhs) const
    { return mult(rhs); }

    T const& operator[] (std::size_t i) const
    { return m[i]; }

    T& operator[] (std::size_t i)
    { return m[i]; }

    Matrix4<T>& operator= (Matrix4<T> const& rhs)
    {
      for (std::size_t i = 0; i < 16; ++i) m[i] = rhs[i];
      return *this;
    }

    //Matrix4<T>& operator, (T const& rhs)
    //{}

    friend std::ostream& operator<< (std::ostream& os, Matrix4<T> const& x)
    {
      return os << "[" << x[0] << "," << x[1] << "," << x[2] << "," << x[3]
          << "; " << x[4] << "," << x[5] << "," << x[6] << "," << x[7]
          << "; " << x[8] << "," << x[9] << "," << x[10] << "," << x[11]
          << "; " << x[12] << "," << x[13] << "," << x[14] << "," << x[15]
          << "]";
    }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_MATRIX4_HEADER */
