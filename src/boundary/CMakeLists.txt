aux_source_directory(. BOUNDARY_SRC)
add_library(xfdtd_boundary OBJECT ${BOUNDARY_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_boundary> PARENT_SCOPE)