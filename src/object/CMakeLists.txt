aux_source_directory(. XFDTD_OBJECT_SRC)
add_library(xfdtd_object OBJECT ${XFDTD_OBJECT_SRC} )
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_object> PARENT_SCOPE)