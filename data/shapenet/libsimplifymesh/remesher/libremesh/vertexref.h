#ifndef REMESHER_VERTEX_REF_HEADER
#define REMESHER_VERTEX_REF_HEADER

#include <vector>

#include "defines.h"
#include "refptrarray.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "vec2.h"

REMESHER_NAMESPACE_BEGIN

struct VertexRef
{
  std::size_t face;
  Vec2f bary;

  VertexRef (void);
  VertexRef (std::size_t face, Vec2f const& bary);
};

/* ---------------------------------------------------------------- */

/*
 * This class is a proxy class for a vector that holds a reference
 * point on a reference mesh for each vertex. The reference point
 * is given with a face index and a barycentric coordinate for that
 * face. This uniquely defines a point on the reference mesh.
 */

class VertexRefList;
typedef RefPtrArray<VertexRefList, VertexRef> VertexRefListPtr;

class VertexRefList : public std::vector<VertexRef>
{
  protected:
    VertexRefList (void);

  public:
    static VertexRefListPtr create (void);
    VertexRefListPtr create_copy (void) const;

    /* Initializes the vertex references with identity mapping. */
    void init_identity (TriangleMeshPtr mesh, VertexInfoListPtr vinfo);

    std::size_t get_memory_usage (void) const;
};

/* ---------------------------------------------------------------- */

inline
VertexRef::VertexRef (void)
{
}

inline
VertexRef::VertexRef (std::size_t face, Vec2f const& bary)
  : face(face), bary(bary)
{
}

inline VertexRefListPtr
VertexRefList::create (void)
{
  return VertexRefListPtr(new VertexRefList);
}

inline
VertexRefList::VertexRefList (void)
{
}

inline VertexRefListPtr
VertexRefList::create_copy (void) const
{
  return VertexRefListPtr(new VertexRefList(*this));
}

inline std::size_t
VertexRefList::get_memory_usage (void) const
{
  return this->capacity() * sizeof(VertexRef);
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_VERTEX_REF_HEADER */
