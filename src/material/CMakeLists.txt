add_library(xfdtd_material OBJECT material.cpp dispersive_material.cpp)
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_material> PARENT_SCOPE)
