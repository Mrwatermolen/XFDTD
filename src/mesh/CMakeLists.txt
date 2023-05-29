aux_source_directory(. MESH_SRC)
add_library(xfdtd_mesh OBJECT ${MESH_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_mesh> PARENT_SCOPE)
