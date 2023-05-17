aux_source_directory(. MONITOR_SRC)
add_library(xfdtd_monitor OBJECT ${MONITOR_SRC})
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_monitor> PARENT_SCOPE)
