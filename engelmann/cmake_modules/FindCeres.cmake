# - Try to find Ceres
# Once done, this will define
#
#  Ceres_FOUND - system has Ceres
#  Ceres_INCLUDE_DIRS - the Ceres include directories
#  Ceres_LIBRARIES - link these to use Ceres

include(LibFindMacros)

# Include dir
find_path(Ceres_INCLUDE_DIR
  NAMES ceres/ceres.h
  PATHS
  /usr/local/include
  PATH_SUFFIXES ceres
)

# Libraries
find_library(Ceres_LIBRARY
  NAMES ceres
  PATHS
  /usr/local/lib
  NO_CMAKE_SYSTEM_PATH
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(Ceres_PROCESS_INCLUDES Ceres_INCLUDE_DIR)
set(Ceres_PROCESS_LIBRARIES Ceres_LIBRARY)
libfind_process(Ceres)
