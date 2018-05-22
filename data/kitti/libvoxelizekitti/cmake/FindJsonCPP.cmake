# - Try to find JsonCPP
# Once done, this will define
#
#  JsonCPP_FOUND - system has JsonCPP
#  JsonCPP_INCLUDE_DIRS - the JsonCPP include directories
#  JsonCPP_LIBRARIES - link these to use JsonCPP

include(LibFindMacros)

# Include dir
find_path(JsonCPP_INCLUDE_DIR
  NAMES json/json.h
  PATHS
  /usr/local/include
  /usr/include
  /usr/local/include
  PATH_SUFFIXES jsoncpp
)

# Libraries
find_library(JsonCPP_LIBRARY
  NAMES jsoncpp
  PATHS
  /usr/lib
  /usr/local/lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this lib depends on.
set(JsonCPP_PROCESS_INCLUDES JsonCPP_INCLUDE_DIR)
set(JsonCPP_PROCESS_LIBRARIES JsonCPP_LIBRARY)
libfind_process(JsonCPP)
