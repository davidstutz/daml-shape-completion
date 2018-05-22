#ifndef REMESHER_ANGLE_SMOOTHING_HEADER
#define REMESHER_ANGLE_SMOOTHING_HEADER

#include "defines.h"
#include "relaxation.h"

REMESHER_NAMESPACE_BEGIN

class AngleSmoothing : public Relaxation
{
  protected:
    Vec2f micropatch_relocate (MicroDelaunay const& md);

  public:
    AngleSmoothing (void);
};

/* ---------------------------------------------------------------- */

inline
AngleSmoothing::AngleSmoothing (void)
{
}

REMESHER_NAMESPACE_END

#endif /* REMESHER_ANGLE_BASED_SMOOTHING_HEADER */
