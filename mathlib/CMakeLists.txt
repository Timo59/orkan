add_library(Math OBJECT linalg.c)

target_include_directories(Math INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set(install_signal OLD OMP PROFILE)
if (${CIRCMETH} IN_LIST install_signal OR ${CIRCMETH} IN_LIST install_signal)
    install(FILES linalg.h DESTINATION include/qlib)
endif ()