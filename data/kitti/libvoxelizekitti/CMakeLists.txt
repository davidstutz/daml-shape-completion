cmake_minimum_required(VERSION 3.2)
project(libvoxelizekitti)

set(CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS} -O3 -g")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(JsonCPP REQUIRED)
find_package(Eigen3 REQUIRED)

message("Eigen3: ${EIGEN3_INCLUDE_DIR}")

include_directories(
  ${JsonCPP_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)
add_executable(voxelize_kitti voxelize_kitti.cpp box/box.cc icp/icp.cpp icp/icpPointToPlane.cpp icp/icpPointToPoint.cpp icp/kdtree.cpp icp/matrix.cpp)
target_link_libraries(voxelize_kitti
  ${JsonCPP_LIBRARIES}
  ${Boost_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
)