aux_source_directory(. XFDTD_EMF_SRC)
add_library(xfdtd_electromagnetic_field OBJECT ${XFDTD_EMF_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_electromagnetic_field> PARENT_SCOPE)