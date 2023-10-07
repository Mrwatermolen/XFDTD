aux_source_directory(. XFDTD_FDTD_BASIC_COFF_SRC)
add_library(xfdtd_fdtd_basic_coff OBJECT ${XFDTD_FDTD_BASIC_COFF_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_fdtd_basic_coff> PARENT_SCOPE)