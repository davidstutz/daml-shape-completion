#ifndef REMESHER_AREA_BASED_RELAXATION
#define REMESHER_AREA_BASED_RELAXATION

#include "defines.h"
#include "vec2.h"
#include "microdelaunay.h"
#include "relaxation.h"

REMESHER_NAMESPACE_BEGIN

class AreaBasedRelaxation : public Relaxation
{
  protected:
    Vec2f micropatch_relocate (MicroDelaunay const& md);

  public:
    AreaBasedRelaxation (void);
};

/* ---------------------------------------------------------------- */

inline
AreaBasedRelaxation::AreaBasedRelaxation (void)
{
}

inline Vec2f
AreaBasedRelaxation::micropatch_relocate (MicroDelaunay const& md)
{
  return md.get_area_optimal_center();
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_AREA_BASED_RELAXATION */
