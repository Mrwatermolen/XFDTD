add_executable(test_1d_fdtd test_1d_fdtd.cpp)
add_executable(test_shape test_shape.cpp)
add_executable(test_source test_source.cpp)
add_executable(test_tensor Test_tensor.cpp)
add_executable(example_fdtd_1d_cpml example_fdtd_1d_cpml.cpp)

target_link_libraries(test_1d_fdtd PRIVATE shape material object waveform
                                      additive_source simulation boundary)
target_link_libraries(test_shape PRIVATE shape)
target_link_libraries(test_source PRIVATE waveform)
target_link_libraries(example_fdtd_1d_cpml PRIVATE shape  waveform additive_source )