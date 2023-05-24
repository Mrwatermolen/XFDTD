aux_source_directory(. EMF_SRC)
add_library(xfdtd_electromagnetic_field OBJECT ${EMF_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_electromagnetic_field> PARENT_SCOPE)