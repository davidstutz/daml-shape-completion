# - Try to find Eigen3
# Once done, this will define
#
#  Eigen3_FOUND - system has Eigen3
#  Eigen3_INCLUDE_DIRS - the Eigen3 include directories
#  Eigen3_LIBRARIES - link these to use Eigen3

include(LibFindMacros)

# Include dir
find_path(Eigen3_INCLUDE_DIR
  NAMES signature_of_eigen3_matrix_library
  PATHS
  /usr/local/include
  PATH_SUFFIXES eigen3
  NO_CMAKE_SYSTEM_PATH
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(Eigen3_PROCESS_INCLUDES Eigen3_INCLUDE_DIR)
libfind_process(Eigen3)
