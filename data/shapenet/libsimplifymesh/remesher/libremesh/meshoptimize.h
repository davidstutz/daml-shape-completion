#ifndef REMESHER_MESH_OPTIMIZE_HEADER
#define REMESHER_MESH_OPTIMIZE_HEADER

#include <vector>
#include <stack>
#include <set>

#include "trianglemesh.h"
#include "vertexinfo.h"
#include "defines.h"

REMESHER_NAMESPACE_BEGIN

/*
 * The algorithm in this class optimizes a mesh by reordering
 * faces and vertices for vertex locality. This makes efficient
 * use of cache, e.g. to accelerate rendering on the graphics card
 * and for the dynamic patch-wise appraoch that uses patch caching.
 * The implementation is based on the method described in the paper:
 *
 *   Pedro V. Sander, Diego Nehab, Joshua Barczak
 *   Fast Triangle Reordering for Vertex Locality and Reduced Overdraw
 *
 * The algorithm is very fast, runs in linear time with linear
 * memory usage and can even be used interactively.
 */

class MeshOptimize
{
  public:
    typedef std::vector<std::size_t> VertexIndexRelocList;

  private:
    typedef std::vector<std::size_t> VertexLiveFaceCounts;
    typedef std::vector<std::size_t> VertexCacheTimes;
    typedef std::stack<std::size_t> DeadEndStack;
    typedef std::vector<bool> FaceEmittedFlags;
    typedef std::set<std::size_t> VertexCandidates;

  private:
    TriangleMeshPtr mesh;
    VertexInfoListPtr vinfo;
    VertexIndexRelocList* irlist;

  private:
    void optimize_face_ordering (std::size_t cache_size);
    void optimize_vertex_ordering (void);

    std::size_t get_next_vertex (std::size_t& cur, std::size_t cache_size,
        VertexCandidates const& next, VertexCacheTimes const& cts,
        std::size_t ts, VertexLiveFaceCounts const& live,
        DeadEndStack& dead);

    std::size_t skip_dead_end (VertexLiveFaceCounts const& live,
        DeadEndStack& dead, std::size_t& cur);

  public:
    MeshOptimize (void);

    /* Sets mesh and vertex info, which is needed to operate. */
    void set_mesh (TriangleMeshPtr mesh);
    void set_vertex_info (VertexInfoListPtr vinfo);

    /* Requests an index relocation list, which gives a permutation
     * of the vertex indices, mapping old positions to new positions. */
    void fill_vertex_reloc_list (VertexIndexRelocList& irlist);
    void unset_vertex_reloc_list (void);

    void optimize (std::size_t cache_size = 24);
    void randomize (void);
};

/* ---------------------------------------------------------------- */

inline
MeshOptimize::MeshOptimize (void)
{
  this->irlist = 0;
}

inline void
MeshOptimize::set_mesh (TriangleMeshPtr mesh)
{
  this->mesh = mesh;
}

inline void
MeshOptimize::set_vertex_info (VertexInfoListPtr vinfo)
{
  this->vinfo = vinfo;
}

inline void
MeshOptimize::fill_vertex_reloc_list (VertexIndexRelocList& irlist)
{
  this->irlist = &irlist;
}

inline void
MeshOptimize::unset_vertex_reloc_list (void)
{
  this->irlist = 0;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_MESH_OPTIMIZE_HEADER */
