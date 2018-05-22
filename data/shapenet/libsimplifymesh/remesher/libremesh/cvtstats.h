#ifndef REMESHER_CVT_STATS_HEADER
#define REMESHER_CVT_STATS_HEADER

#include "trianglemesh.h"
#include "densityfield.h"
#include "vertexref.h"
#include "vertexinfo.h"
#include "defines.h"

REMESHER_NAMESPACE_BEGIN

class CvtStats
{
  private:
    TriangleMeshPtr rmesh;
    DensityFieldPtr rdens;

    TriangleMeshPtr emesh;
    VertexInfoListPtr evinfo;
    VertexRefListPtr evref;

  public:
    void set_ref_mesh (TriangleMeshPtr rmesh);
    void set_ref_density (DensityFieldPtr rdens);

    void set_evo_mesh (TriangleMeshPtr emesh);
    void set_evo_vinfo (VertexInfoListPtr evinfo);
    void set_evo_vref (VertexRefListPtr evref);

    void calc_stats (std::string const& gp_outfile);
};

/* ---------------------------------------------------------------- */

inline void
CvtStats::set_ref_mesh (TriangleMeshPtr rmesh)
{
  this->rmesh = rmesh;
}

inline void
CvtStats::set_ref_density (DensityFieldPtr rdens)
{
  this->rdens = rdens;
}

inline void
CvtStats::set_evo_mesh (TriangleMeshPtr emesh)
{
  this->emesh = emesh;
}

inline void
CvtStats::set_evo_vinfo (VertexInfoListPtr evinfo)
{
  this->evinfo = evinfo;
}

inline void
CvtStats::set_evo_vref (VertexRefListPtr evref)
{
  this->evref = evref;
}


REMESHER_NAMESPACE_END

#endif /* REMESHER_CVT_STATS_HEADER */
