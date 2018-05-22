#ifndef REMESHER_MODEL_LOADER_HEADER
#define REMESHER_MODEL_LOADER_HEADER

#include <string>

#include "defines.h"
#include "trianglemesh.h"

REMESHER_NAMESPACE_BEGIN

class ModelLoader
{
  public:
    static TriangleMeshPtr load_model (std::string const& filename);

    static TriangleMeshPtr load_off_model (std::string const& filename);
    static TriangleMeshPtr load_ply_model (std::string const& filename);
    static TriangleMeshPtr load_obj_model (std::string const& filename);
    static TriangleMeshPtr load_m_model (std::string const& filename);
};

REMESHER_NAMESPACE_END

#endif /* REMESHER_MODEL_LOADER_HEADER */
