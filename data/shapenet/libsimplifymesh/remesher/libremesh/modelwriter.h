#ifndef REMESHER_MODEL_WRITER_HEADER
#define REMESHER_MODEL_WRITER_HEADER

#include <string>

#include "defines.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

class ModelWriter
{
  public:
    /* Save model by filename extension. */
    static void save_model (std::string const& filename,
        TriangleMeshPtr mesh, bool with_normals = false);

    /* Save model regardless of filename. */
    static void save_off_model (std::string const& filename,
        TriangleMeshPtr mesh, bool with_normals = false);
    static void save_ply_model (std::string const& filename,
        TriangleMeshPtr mesh, bool with_normals = false);
    static void save_svg_file (std::string const& filename,
        TriangleMeshPtr mesh);
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_MODEL_WRITER_HEADER */
