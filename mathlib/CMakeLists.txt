add_library(Math OBJECT linalg.c)

target_include_directories(Math INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(FILES linalg.h DESTINATION include)