# src/CMakeLists.txt - Main library target

# Collect source files
set(Q_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/state.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/state_blocked.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/gate.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/qhipster.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/mhipster.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/mhipster_block.c"
)

# Include directory
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# Create library target
add_library(q SHARED ${Q_SOURCES})

set_target_properties(q PROPERTIES
        VERSION         ${PROJECT_VERSION}
        SOVERSION       ${PROJECT_VERSION_MAJOR}
        PUBLIC_HEADER   "${INCLUDE_DIR}/qlib.h"
)
target_include_directories(q
        PUBLIC
            $<BUILD_INTERFACE:${INCLUDE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(q
        PUBLIC
            ${PROJECT_COMPILER_FLAGS}
            blas_compiler_flags
            omp_compiler_flags
)