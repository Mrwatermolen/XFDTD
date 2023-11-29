cmake_minimum_required(VERSION 3.20)

project(XFDTD_EXAMPLE VERSION 0.0.0 LANGUAGES CXX)

set(XFDTD_EXAMPLE_MAIN_PROJECT OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(XFDTD_EXAMPLE_MAIN_PROJECT ON)
endif()

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(THIRD_PARTY_LIB_XFDTD xfdtd)

add_executable(xfdtd_example_2d_cylinde_scatter 2d_cylinder_scatter.cpp)
target_link_libraries(xfdtd_example_2d_cylinde_scatter ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_3d_pec_sphere_monostatic_rcs 3d_pec_sphere_monostatic_rcs.cpp)
target_link_libraries(xfdtd_example_3d_pec_sphere_monostatic_rcs ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_3d_pec_sphere_bistatic_rcs 3d_pec_sphere_bistatic_rcs.cpp)
target_link_libraries(xfdtd_example_3d_pec_sphere_bistatic_rcs ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_basic_lumped_element basic_lumped_element.cpp)
target_link_libraries(xfdtd_example_basic_lumped_element ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_capacitor capacitor.cpp)
target_link_libraries(xfdtd_example_capacitor ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_inductor inductor.cpp)
target_link_libraries(xfdtd_example_inductor ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_rlc_circuit rlc_circuit.cpp)
target_link_libraries(xfdtd_example_rlc_circuit ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_s_parameter_network s_parameter_network.cpp)
target_link_libraries(xfdtd_example_s_parameter_network ${THIRD_PARTY_LIB_XFDTD})

add_executable(xfdtd_example_thin_wire thin_wire.cpp)
target_link_libraries(xfdtd_example_thin_wire ${THIRD_PARTY_LIB_XFDTD})
