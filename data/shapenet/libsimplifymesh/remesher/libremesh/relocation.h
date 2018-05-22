#ifndef REMESHER_RELOCATION_HEADER
#define REMESHER_RELOCATION_HEADER

#include <list>

#include "defines.h"
#include "trianglemesh.h"
#include "vertexinfo.h"
#include "vertexref.h"
#include "triangle2.h"
#include "patchcache.h"

#if REMESHER_FB_PARAMETRIZATION
#   include "patch2d_fb.h"
#else
#   include "patch2d.h"
#endif

REMESHER_NAMESPACE_BEGIN

/*
 * These classes are capable of relocating a single vertex. Given three
 * point references and a barycentric coordinate regaring the three
 * points, the algorithm finds a triangle on the reference mesh that
 * contains the new point and a barycentric coordinate that uniquely
 * defines the new point. The class implements a technique from:
 *
 *   Vitaly Surazhsky and Craig Gotsman
 *   Explicit Surface Remeshing
 *
 * The actual relocation is done using a parametrization created
 * with the Patch2d class.
 */

/* ---------------------------------------------------------------- */

class Relocation
{
#if REMESHER_PATCH_CACHING
  private:
    PatchCachePtr cache;
#endif

  private:
    /* Reference mesh data structures. */
    TriangleMeshPtr rmesh;
    VertexInfoListPtr vinfo;

  private:
    Patch2dPtr get_patch (std::size_t f1, std::size_t f2, std::size_t f3);
    std::size_t find_opposite_face (std::size_t face, std::size_t fvidx);

  public:
    Relocation (void);
    ~Relocation (void);

    void set_data (TriangleMeshPtr ref_mesh, VertexInfoListPtr ref_vinfo);

    VertexRef relocate (VertexRef const& v1, VertexRef const& v2,
        VertexRef const& v3, Vec2f const& bary_coords);
};

/* ---------------------------------------------------------------- */

inline
Relocation::Relocation (void)
{
#if REMESHER_PATCH_CACHING
  this->cache = PatchCache::create(REMESHER_PATCH_CACHE_SIZE);
#endif
}

inline
Relocation::~Relocation (void)
{
  this->cache->print_statistics();
}

inline void
Relocation::set_data (TriangleMeshPtr ref_mesh, VertexInfoListPtr ref_vinfo)
{
  this->rmesh = ref_mesh;
  this->vinfo = ref_vinfo;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_RELOCATION_HEADER */
