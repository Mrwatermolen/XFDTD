add_library(object OBJECT object.cpp)

target_link_libraries(object
   PRIVATE 
   shape   material )