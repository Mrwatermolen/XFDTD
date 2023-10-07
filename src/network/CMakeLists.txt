aux_source_directory(. XFDTD_NETWORK_SRC)
add_library(xfdtd_network OBJECT ${XFDTD_NETWORK_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_network> PARENT_SCOPE)