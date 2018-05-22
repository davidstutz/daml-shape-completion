#include <iostream>

#include "patchcache.h"

REMESHER_NAMESPACE_BEGIN

PatchCache::PatchCache (void)
{
  this->cache_size = 0;
  this->max_cache_size = REMESHER_PATCH_CACHE_SIZE;

#if REMESHER_FAST_PATCH_LOOKUP
  this->clock = 0;
#endif

  this->clear_statistics();
}

/* ---------------------------------------------------------------- */

PatchCache::~PatchCache (void)
{
  // std::cout << "Destroying cache with " << (this->cache_size / 1024)
      // << " KB of cached patches" << std::endl;

#if REMESHER_FAST_PATCH_LOOKUP
  for (std::size_t i = 0; i < this->patches.size(); ++i)
    delete this->patches[i];
#endif
}

/* ---------------------------------------------------------------- */

Patch2dPtr
PatchCache::lookup_patch (std::size_t f1, std::size_t f2, std::size_t f3)
{
  /* Record cache hit/miss statistics. Assume a hit and remove at miss. */
  this->cache_hits += 1;

#if REMESHER_PARALLELIZATION
  /* Aquire a read lock on the cache if parallelization is enabled. */
  ReadLock rlock(this->cache_lock);
#endif

#if REMESHER_FAST_PATCH_LOOKUP

  /*
   * Use fast but memory intensive patch lookup on the data structure.
   * We select a patch that is present for all three faces.
   */
  std::size_t cs = this->cache.size();
  if (cs > f1 && cs > f2 && cs > f3)
  {
    CachedPatches const& p1 = this->cache[f1];
    CachedPatches const& p2 = this->cache[f2];
    CachedPatches const& p3 = this->cache[f3];

    for (std::size_t i = 0; i < p1.size(); ++i)
    {
      if (p1[i] == 0)
        continue;

      for (std::size_t j = 0; j < p2.size(); ++j)
      {
        if (p1[i] != p2[j])
          continue;

        for (std::size_t k = 0; k < p3.size(); ++k)
          if (p1[i] == p3[k])
          {
            this->clock += 1;
            p1[i]->second = this->clock;
            return p1[i]->first;
          }
      }
    }
  }

#else /* No fast cache lookup. */

  /*
   * Use a slower but memory-free patch lookup on the data structure.
   * This technique iterates over all cached patches and tries to find
   * a patch that contains all three faces. This technique is not suited
   * for parallelization, because cache hits should be moved to the front.
   */
  PatchCacheDS::iterator iter;
  for (iter = this->cache.begin(); iter != this->cache.end(); iter++)
  {
    std::size_t f1_2d = (*iter)->lookup_face(f1);
    if (f1_2d == MAX_SIZE_T)
      continue;

    std::size_t f2_2d = (*iter)->lookup_face(f2);
    if (f2_2d == MAX_SIZE_T)
      continue;

    std::size_t f3_2d = (*iter)->lookup_face(f3);
    if (f3_2d == MAX_SIZE_T)
      continue;

    Patch2dPtr patch2d = *iter;

# if !REMESHER_PARALLELIZATION
    /* If we have a cache hit, we improve caching by moving the patch to
     * the front in the cache data structure. This, however, is not possible
     * with parallelization, because every successful lookup would also
     * require a write operation on the cache (a write look). */
    if (iter != this->cache.begin())
    {
      this->cache.erase(iter);
      this->cache.push_front(patch2d);
    }
# endif

    return patch2d;
  }
#endif

  /* Record statistics. Remove assumed hit and count the miss. */
  this->cache_hits -= 1;
  this->cache_misses += 1;

  /* Return an empty patch to indicate a cache miss. */
  return Patch2dPtr();
}

/* ---------------------------------------------------------------- */

void
PatchCache::cache_patch (Patch2dPtr patch2d)
{
# if REMESHER_PARALLELIZATION
  //MutexLock lock(this->cache_lock);
  WriteLock wlock(this->cache_lock);
# endif

  /* Get the size of the patch. */
  std::size_t patch_size = patch2d->get_memory_usage();

  /* Record statistics. Keeping track of deletions
   * is up to the specific algorithms below. */
  this->p_created += 1;
  this->p_created_size += patch_size;


#if REMESHER_FAST_PATCH_LOOKUP

  /* Remove old patches if the cache size is going to be exceeded. */
  if (this->cache_size + patch_size > this->max_cache_size
      && !this->cache.empty())
  {
    /* Remove 20% of the oldest patches (factor 0.8). */
    std::size_t del_idx = (std::size_t)((float)this->patches.size() * 0.8f);
    // std::cout << "Cache full, removing patches, "
        // << this->patches.size() << " - "
        // << (this->patches.size() - del_idx)
        // << " = " << del_idx << std::endl;

    /* We first sort the patches by last access. */
    std::sort(this->patches.begin(), this->patches.end(),
        PatchCache::entry_cmp);

    /* Remove the patches from the DS and shrink the vector. */
    for (std::size_t i = del_idx; i < this->patches.size(); ++i)
    {
      CacheEntry* to_delete = this->patches[i];
      Patch2d::FaceMapping const& fmap = to_delete->first->get_face_mapping();
      for (std::size_t j = 0; j < fmap.size(); ++j)
      {
        CachedPatches& cp = this->cache[fmap[j]];
        for (std::size_t k = 0; k < cp.size(); ++k)
          if (cp[k] == to_delete)
            cp[k] = 0;
      }
      std::size_t free_size = to_delete->first->get_memory_usage();
      this->cache_size -= free_size;
      delete to_delete;

      /* Record statistics. */
      this->p_deleted += 1;
      this->p_deleted_size += free_size;
    }
    this->patches.resize(del_idx);
  }

  /* Appropriately resize the face cache. We use the last face index
   * because it's the highest index in the face mapping. */
  Patch2d::FaceMapping const& fmap = patch2d->get_face_mapping();
  if (!fmap.empty())
    if (fmap.back() >= this->cache.size())
      this->cache.resize(fmap.back() + 1);

  /* Add the new patch to the cache. */
  this->clock += 1;
  CacheEntry* new_entry = new CacheEntry;
  new_entry->first = patch2d;
  new_entry->second = this->clock;
  this->patches.push_back(new_entry);
  for (std::size_t i = 0; i < fmap.size(); ++i)
  {
    CachedPatches& cp = this->cache[fmap[i]];
    bool inserted = false;
    for (std::size_t j = 0; j < cp.size(); ++j)
      if (cp[j] == 0)
      {
        cp[j] = new_entry;
        inserted = true;
      }

    if (!inserted)
      cp.insert(cp.begin(), new_entry);
  }

#else /* No fast cache lookup. */

  /* Write 2D patch to the cache. */
  while (this->cache_size + patch_size > this->max_cache_size
      && !this->cache.empty())
  {
    //std::cout << "Cache full, removing patch" << std::endl;
    std::size_t free_size = this->cache.back()->get_memory_usage();
    this->cache.pop_back();
    this->cache_size -= free_size;

    /* Record statistics. */
    this->p_deleted += 1;
    this->p_deleted_size += free_size;
  }
  this->cache.push_front(patch2d);

#endif

  /* Add the patch size to the current cache size. */
  this->cache_size += patch_size;
}

/* ---------------------------------------------------------------- */

void
PatchCache::print_statistics (void)
{
  std::size_t cache_kb = this->cache_size / 1024;
  std::size_t max_cache_kb = this->max_cache_size / 1024;

  // std::cout << "Patch caching statistics" << std::endl
      // << "  Current usage: " << cache_kb << " KB of " << max_cache_kb
      // << " KB (" << 100 * cache_kb / max_cache_kb << "%)" << std::endl
      // << "  Patches created: " << this->p_created
      // << " (" << this->p_created_size / 1024 << " KB)" << std::endl
      // << "  Patches deleted: " << this->p_deleted
      // << " (" << this->p_deleted_size / 1024 << " KB)" << std::endl
      // << "  Cache hits: " << this->cache_hits
      // << ", Cache misses: " << this->cache_misses
      // << ", Ratio: " << ((float)this->cache_hits / (float)this->cache_misses)
      // << std::endl;

  std::size_t overhead = 0;

#if REMESHER_FAST_PATCH_LOOKUP
  overhead += this->cache.capacity() * sizeof(CachedPatches);
  std::size_t patch_refs = 0;
  for (std::size_t i = 0; i < this->cache.size(); ++i)
  {
    overhead += this->cache[i].capacity() * sizeof(CacheEntry*);
    patch_refs += this->cache[i].size();
  }
  overhead += this->patches.capacity() * sizeof(CacheEntry*);
  overhead += this->patches.size() * sizeof(CacheEntry);
  // std::cout << "  Average patches per triangle: "
      // << (float)patch_refs / (float)this->cache.size() << std::endl;
#else
  overhead += this->cache.size() * sizeof(Patch2dPtr);
#endif

  // std::cout << "  DS overhead: " << (overhead / 1024) << " KB" << std::endl;
}

/* ---------------------------------------------------------------- */

void
PatchCache::clear_statistics (void)
{
  this->p_created = 0;
  this->p_created_size = 0;
  this->p_deleted = 0;
  this->p_deleted_size = 0;
  this->cache_hits = 0;
  this->cache_misses = 0;
}

REMESHER_NAMESPACE_END
