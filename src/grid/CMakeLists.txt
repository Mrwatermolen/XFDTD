aux_source_directory(. XFDTD_GRID_SRC)
add_library(xfdtd_grid OBJECT ${XFDTD_GRID_SRC})
set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_grid>
    PARENT_SCOPE)
