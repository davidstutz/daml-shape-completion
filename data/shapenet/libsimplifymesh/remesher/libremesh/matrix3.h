#ifndef REMESHER_MATRIX3_HEADER
#define REMESHER_MATRIX3_HEADER

#include <ostream>

#include "defines.h"
#include "vec3.h"

REMESHER_NAMESPACE_BEGIN

template <class T> class Matrix3;

typedef Matrix3<float> Matrix3f;
typedef Matrix3<double> Matrix3d;
typedef Matrix3<int> Matrix3i;
typedef Matrix3<unsigned int> Matrix3ui;
typedef Matrix3<char> Matrix3c;

template <class T>
class Matrix3
{
  public:
    T m[9];

    Matrix3 (void) {}
    Matrix3 (T* matrix)
    { for (std::size_t i = 0; i < 9; ++i) m[i] = matrix[i]; }

    Matrix3<T> plus (Matrix3<T> const& rhs) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i) ret[i] = m[i] + rhs[i];
      return ret;
    }

    Matrix3<T> minus (Matrix3<T> const& rhs) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i) ret[i] = m[i] - rhs[i];
      return ret;
    }

    Matrix3<T> mult (Matrix3<T> const& rhs) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i)
      {
        ret[i] = (T)0;
        for (std::size_t j = 0; j < 3; ++j)
          ret[i] += m[i-i%3 + j] * rhs[i%3 + j*3];
      }
      return ret;
    }

    Matrix3<T> mult (T const& rhs) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i) ret[i] = m[i] * rhs;
      return ret;
    }

    Vec3<T> mult (Vec3<T> const& rhs) const
    {
      return Vec3<T>(m[0] * rhs.x + m[1] * rhs.y + m[2] * rhs.z,
          m[3] * rhs.x + m[4] * rhs.y + m[5] * rhs.z,
          m[6] * rhs.x + m[7] * rhs.y + m[8] * rhs.z);
    }

    Matrix3<T> div (T const& rhs) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i) ret[i] = m[i] / rhs;
      return ret;
    }

    Matrix3<T> neg (void) const
    {
      Matrix3<T> ret;
      for (std::size_t i = 0; i < 9; ++i) ret[i] = -m[i];
      return ret;
    }

    Matrix3<T> transpose (void) const
    {
      Matrix3<T> ret;
      ret[0] = m[0]; ret[1] = m[3]; ret[2] = m[6];
      ret[3] = m[1]; ret[4] = m[4]; ret[5] = m[7];
      ret[6] = m[2]; ret[7] = m[5]; ret[8] = m[8];
      return ret;
    }

    T const* array (void) const
    { return m; }

    T* array (void)
    { return m; }

    Matrix3<T> operator+ (Matrix3<T> const& rhs) const
    { return plus(rhs); }

    Matrix3<T> operator- (Matrix3<T> const& rhs) const
    { return minus(rhs); }

    Matrix3<T> operator* (Matrix3<T> const& rhs) const
    { return mult(rhs); }

    Matrix3<T> operator* (T const& rhs) const
    { return mult(rhs); }

    Vec3<T> operator* (Vec3<T> const& rhs) const
    { return mult(rhs); }

    T const& operator[] (std::size_t i) const
    { return m[i]; }

    T& operator[] (std::size_t i)
    { return m[i]; }

    Matrix3<T>& operator= (Matrix3<T> const& rhs)
    {
      for (std::size_t i = 0; i < 9; ++i) m[i] = rhs[i];
      return *this;
    }

    //Matrix3<T>& operator, (T const& rhs)
    //{}

    friend std::ostream& operator<< (std::ostream& os, Matrix3<T> const& x)
    {
      return os << "[" << x[0] << "," << x[1] << "," << x[2]
          << "; " << x[3] << "," << x[4] << "," << x[5]
          << "; " << x[6] << "," << x[7] << "," << x[8]
          << "]";
    }
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_MATRIX3_HEADER */
