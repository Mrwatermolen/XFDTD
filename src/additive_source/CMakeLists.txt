aux_source_directory(. DIR_ADDITIVE_SOURCE_SRCS)
add_library(additive_source OBJECT ${DIR_ADDITIVE_SOURCE_SRCS})

target_link_libraries(additive_source PRIVATE waveform)
