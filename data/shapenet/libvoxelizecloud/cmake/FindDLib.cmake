# - Try to find DLIB
# Once done, this will define
#
#  DLIB_FOUND - system has DLIB
#  DLIB_INCLUDE_DIRS - the DLIB include directories
#  DLIB_LIBRARIES - link these to use DLIB

include(LibFindMacros)

# Include dir
find_path(DLIB_INCLUDE_DIR
  NAMES dlib/config.h
  PATHS
  /usr/local/include
  PATH_SUFFIXES dlib
  # https://stackoverflow.com/questions/16875664/preventing-cmake-from-finding-installed-libraries-instead-of-local-libraries
  NO_CMAKE_SYSTEM_PATH
)

# Libraries
find_library(DLIB_LIBRARY
  NAMES dlib
  PATHS
  /usr/lib
  /usr/local/lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(DLIB_PROCESS_INCLUDES DLIB_INCLUDE_DIR)
set(DLIB_PROCESS_LIBRARIES DLIB_LIBRARY)
libfind_process(DLIB)
