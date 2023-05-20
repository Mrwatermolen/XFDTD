aux_source_directory(. TFSF_SRC)
add_library(xfdtd_tfsf OBJECT ${TFSF_SRC})
set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_tfsf>
    PARENT_SCOPE)
