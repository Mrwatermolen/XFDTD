aux_source_directory(. XFDTD_BOUNDARY_SRC)
add_library(xfdtd_boundary OBJECT ${XFDTD_BOUNDARY_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_boundary> PARENT_SCOPE)