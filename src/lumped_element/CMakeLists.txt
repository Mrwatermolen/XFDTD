aux_source_directory(. XFDTD_LUMPED_ELEMENT_SRCS)
add_library(xfdtd_lumped_element ${XFDTD_LUMPED_ELEMENT_SRCS})
set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_lumped_element>
    PARENT_SCOPE)