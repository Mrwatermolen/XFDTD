cmake_minimum_required(VERSION 3.20.0)

project(
  XFDTD
  VERSION 0.0.0
  LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed")
endif()

include("cmake/FindFFTW.cmake")
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
set(FFTW_USE_DOUBLE on)
find_package(xtensor-fftw REQUIRED)
if(WIN32)
  find_package(OpenBLAS REQUIRED)
  set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

set(XFDTD_XTENSOR_LIB_INCLUDE_DIR ${xtensor_INCLUDE_DIRS})

include_directories(
  src/include ${FFTW_INCLUDE_DIR} ${XFDTD_XTENSOR_LIB_INCLUDE_DIR}
  ${xtensor-fftw_INCLUDE_DIRS})

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(examples)
