add_library(xfdtd_shape OBJECT shape.cpp sphere.cpp cube.cpp line.cpp)
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_shape> PARENT_SCOPE)
