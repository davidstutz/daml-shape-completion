#ifndef REMESHER_PERMUTE_HEADER
#define REMESHER_PERMUTE_HEADER

#include <vector>

#include "defines.h"
#include "exception.h"

REMESHER_NAMESPACE_BEGIN

/*
 * This class permutes a vector of elements v using a permutation
 * given by vector p, calculating v' = p(v). Vector p can be specified
 * by using a mapping from old indices to new indices,
 *
 *   v'_p[i] = v_i  e.g.  v = [a, b, c], p = [1, 2, 0], v' = [c, a, b].
 *
 * which is better called a index-based relocation of the elements.
 * Each element is copied two times, from the original vector to
 * a temporary variable, and back to the vector.
 * This is available using the "permute_reloc" function.
 *
 * Another way of defining p is mathematically more natural and
 * also computationally more efficient:
 *
 *   v'_i = v_p[i]  e.g.  v = [a, b, c], p = [1, 2, 0], v' = [b, c, a]
 *
 * This is available using the "permute_math" function. Each element
 * is copied only once inside the vector, and elements at the beginning
 * of a cycle are copied twice.
 *
 * "permute_reloc" performs the inverse permutation of "permute_math".
 */

template <class V, class P>
class Permute
{
  public:
    typedef std::vector<V> VecV;
    typedef std::vector<P> VecP;
    typedef std::vector<bool> VecB;

  public:
    static void permute_reloc (VecV& v, VecP const& p);
    static void permute_math (VecV& v, VecP const& p);
};

/* ---------------------------------------------------------------- */

template <class V, class P>
inline void
Permute<V, P>::permute_reloc (typename Permute::VecV& v,
    typename Permute::VecP const& p)
{
  if (v.size() != p.size())
    throw Exception("Invalid permutation arguments");

  if (v.empty())
    return;

  VecB visited(v.size(), false);
  std::size_t i = 0;
  std::size_t seek = 0;
  do
  {
    /* Permute a cycle using index-based relocation permutation. */
    V tmp[2];
    bool idx = false;
    tmp[idx] = v[i];
    while (!visited[i])
    {
      tmp[!idx] = v[p[i]];
      v[p[i]] = tmp[idx];
      idx = !idx;
      visited[i] = true;
      i = p[i];
    }

    /* Seek for an unvisited element. */
    i = MAX_SIZE_T;
    for (; seek < visited.size(); ++seek)
      if (!visited[seek])
      {
        i = seek;
        break;
      }
  }
  while (i != MAX_SIZE_T);
}

/* ---------------------------------------------------------------- */

template <class V, class P>
inline void
Permute<V, P>::permute_math (typename Permute::VecV& v,
    typename Permute::VecP const& p)
{
  if (v.size() != p.size())
    throw Exception("Invalid permutation arguments");

  if (v.empty())
    return;

  VecB visited(v.size(), false);
  std::size_t i = 0;
  std::size_t seek = 0;
  do
  {
    /* Permute a cycle using mathematical permutation. */
    visited[i] = true;
    if (i != p[i])
    {
      V tmp = v[i];
      while (!visited[p[i]])
      {
        v[i] = v[p[i]];
        visited[p[i]] = true;
        i = p[i];
      }
      v[i] = tmp;
    }

    /* Seek for an unvisited element. */
    i = MAX_SIZE_T;
    for (; seek < visited.size(); ++seek)
      if (!visited[seek])
      {
        i = seek;
        break;
      }
  }
  while (i != MAX_SIZE_T);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PERMUTE_HEADER */
