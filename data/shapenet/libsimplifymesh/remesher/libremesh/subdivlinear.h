#ifndef REMEHSER_SUBDIV_LINEAR_HEADER
#define REMEHSER_SUBDIV_LINEAR_HEADER

#include "defines.h"
#include "subdivbase.h"

REMESHER_NAMESPACE_BEGIN

class SubdivLinear : public SubdivBase
{
  protected:
    void subdiv_impl (void);

  public:
    SubdivLinear (void);
};

/* ---------------------------------------------------------------- */

inline
SubdivLinear::SubdivLinear (void)
{
}

REMESHER_NAMESPACE_END

#endif /* REMEHSER_SUBDIV_LINEAR_HEADER */
