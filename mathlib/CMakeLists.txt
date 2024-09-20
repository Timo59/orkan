add_library(Math STATIC linalg.c)

target_include_directories(Math INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)