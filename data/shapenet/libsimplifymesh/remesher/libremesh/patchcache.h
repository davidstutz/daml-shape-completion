#ifndef REMESHER_PATCH_CACHE
#define REMESHER_PATCH_CACHE

#include <list>
#include <vector>

#include "defines.h"
#include "refptr.h"
#include "thread.h"

#if REMESHER_FB_PARAMETRIZATION
#   include "patch2d_fb.h"
#else
#   include "patch2d.h"
#endif

REMESHER_NAMESPACE_BEGIN

class PatchCache;
typedef RefPtr<PatchCache> PatchCachePtr;

class PatchCache
{
  private:
#if REMESHER_FAST_PATCH_LOOKUP
    typedef std::pair<Patch2dPtr, std::size_t> CacheEntry;
    typedef std::vector<CacheEntry*> CachedPatches;
    typedef std::vector<CachedPatches> PatchCacheDS;

    static bool entry_cmp (CacheEntry const* c1, CacheEntry const* c2)
    { return c1->second > c2->second; }

    PatchCacheDS cache;
    CachedPatches patches;
    std::size_t clock;
#else
    typedef std::list<Patch2dPtr> PatchCacheDS;
    PatchCacheDS cache;
#endif

    std::size_t cache_size;
    std::size_t max_cache_size;

#if REMESHER_PARALLELIZATION
    /* The data structure is locked for some cache operations. */
    ReadWriteLock cache_lock;
#endif

    /* Caching statistics. */
    std::size_t p_created;
    std::size_t p_created_size;
    std::size_t p_deleted;
    std::size_t p_deleted_size;
    std::size_t cache_hits;
    std::size_t cache_misses;

  protected:
    PatchCache (void);

  public:
    ~PatchCache (void);
    static PatchCachePtr create (std::size_t max_cache_size);

    Patch2dPtr lookup_patch (std::size_t f1, std::size_t f2, std::size_t f3);
    void cache_patch (Patch2dPtr patch2d);

    void print_statistics (void);
    void clear_statistics (void);
};

/* ---------------------------------------------------------------- */

inline PatchCachePtr
PatchCache::create (std::size_t max_cache_size)
{
  PatchCachePtr ret = PatchCachePtr(new PatchCache);
  ret->max_cache_size = max_cache_size;
  return ret;
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_PATCH_CACHE */
