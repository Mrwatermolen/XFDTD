cmake_minimum_required(VERSION 3.24.0)

project(
  XFdtd
  VERSION 0.0.0
  LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)

set(XFDTD_EIGEN_LIB_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR})

include_directories(src/include ${XFDTD_EIGEN_LIB_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(test)
