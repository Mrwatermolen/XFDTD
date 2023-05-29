aux_source_directory(. NFFFT_SRC)
add_library(xfdtd_nffft OBJECT ${NFFFT_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_nffft> PARENT_SCOPE)
