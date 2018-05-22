#ifndef REMEHSER_SUBDIV_LOOP_HEADER
#define REMEHSER_SUBDIV_LOOP_HEADER

#include "defines.h"
#include "featureedges.h"
#include "subdivbase.h"

REMESHER_NAMESPACE_BEGIN

/*
 * Implements the Loop Subdivision scheme.
 */
class SubdivLoop : public SubdivBase
{
  protected:
    enum LoopVClass
    {
      LVC_SMOOTH,
      LVC_DART,
      LVC_CREASE_REG,
      LVC_CREASE_NONREG,
      LVC_CORNER
    };

  protected:
    std::vector<LoopVClass> vclasses;

  protected:
    void subdiv_impl (void);

    /* Classification of all vertices into LoopVClass classes. */
    void classify_vertices (void);

    /* Three edge subdivision masks. */
    void smooth_subdiv (std::size_t v0idx, std::size_t v1idx,
        std::size_t v2idx, std::size_t v3idx);
    void regular_crease_subdiv(std::size_t v1idx, std::size_t v2idx);
    void nonregular_crease_subdiv(std::size_t v1idx, std::size_t v2idx);

  public:
    SubdivLoop (void);
};

/* ---------------------------------------------------------------- */

inline
SubdivLoop::SubdivLoop (void)
{
}

REMESHER_NAMESPACE_END

#endif /* REMEHSER_SUBDIV_LOOP_HEADER */
