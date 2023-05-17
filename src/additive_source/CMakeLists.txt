aux_source_directory(. DIR_ADDITIVE_SOURCE_SRCS)
add_library(xfdtd_additive_source OBJECT ${DIR_ADDITIVE_SOURCE_SRCS})

set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_additive_source> PARENT_SCOPE)
