aux_source_directory(. DIR_WAVEFORM_SRCS)
add_library(xfdtd_waveform OBJECT ${DIR_WAVEFORM_SRCS})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_waveform> PARENT_SCOPE)