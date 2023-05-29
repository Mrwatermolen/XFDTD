add_subdirectory(boundary)
add_subdirectory(object)
add_subdirectory(material)
add_subdirectory(shape)
add_subdirectory(simulation)
add_subdirectory(waveform)
add_subdirectory(additive_source)
add_subdirectory(monitor)
add_subdirectory(tfsf)
add_subdirectory(electromagnetic_field)
add_subdirectory(nffft)
add_subdirectory(mesh)

add_library(xfdtd STATIC ${ALL_OBJECT_FILES})

set(XFDTD_LIBS
    xfdtd_additive_source
    xfdtd_boundary
    xfdtd_electromagnetic_field
    xfdtd_material
    xfdtd_monitor
    xfdtd_object
    xfdtd_shape
    xfdtd_simulation
    xfdtd_waveform
    xfdtd_nffft)

target_link_libraries(xfdtd ${XFDTD_LIBS} xtensor ${FFTW_LIBRARIES})

target_include_directories(
  xfdtd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)