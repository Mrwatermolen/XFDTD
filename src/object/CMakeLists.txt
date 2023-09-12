add_subdirectory(thin_wire)
add_library(xfdtd_object OBJECT object.cpp)
set(ALL_OBJECT_FILES ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:xfdtd_object> PARENT_SCOPE)