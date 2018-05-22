#ifndef REMESHER_LLOYD_RELAXATION
#define REMESHER_LLOYD_RELAXATION

#include "defines.h"
#include "vec2.h"
#include "microdelaunay.h"
#include "relaxation.h"

REMESHER_NAMESPACE_BEGIN

class LloydRelaxation : public Relaxation
{
  protected:
    Vec2f micropatch_relocate (MicroDelaunay const& md);

  public:
    LloydRelaxation (void);
};

/* ---------------------------------------------------------------- */

inline
LloydRelaxation::LloydRelaxation (void)
{
}

inline Vec2f
LloydRelaxation::micropatch_relocate (MicroDelaunay const& md)
{
  return md.get_voronoi_centroid();
}

REMESHER_NAMESPACE_END

#endif
