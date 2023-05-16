aux_source_directory(. SIMULATION_SRC)
add_library(simulation OBJECT ${SIMULATION_SRC})

target_link_libraries(simulation PRIVATE shape material object waveform)
