cmake_minimum_required(VERSION 2.8)
project(mesh_voxelization)

set(CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
find_package(Eigen3 REQUIRED)

message("Eigen3: ${EIGEN3_INCLUDE_DIR}")

include_directories(${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS})
add_executable(read_hdf5 read_hdf5.cpp)
target_link_libraries(read_hdf5 ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES})

add_executable(read_off read_off.cpp)
target_link_libraries(read_off ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES})

add_executable(read_txt read_txt.cpp)
target_link_libraries(read_txt ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES})