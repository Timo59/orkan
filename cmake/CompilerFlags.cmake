# CompilerFlags.cmake - Compiler and linking configuration

# Set name for interface library holding compiler flags for the project
set(PROJECT_COMPILER_FLAGS "${PROJECT_NAME}_compiler_flags")

# Interface library for project compiler flags
# Optimization and debug flags come from CMAKE_BUILD_TYPE:
#   Debug:   -O0 -g          (cmake default)
#   Release: -O3 -DNDEBUG    (cmake default)
add_library(${PROJECT_COMPILER_FLAGS} INTERFACE)

target_compile_options(${PROJECT_COMPILER_FLAGS}
        INTERFACE
            -fno-strict-aliasing
            $<$<CONFIG:Release>:-march=native>
)
target_link_libraries(${PROJECT_COMPILER_FLAGS}
        INTERFACE $<$<AND:$<BOOL:${UNIX}>,$<NOT:$<BOOL:${APPLE}>>>:-lm>
)

# Interface library to specify usage of BLAS routines
add_library(blas_compiler_flags INTERFACE)

# Use BLAS configuration from Dependencies.cmake
if(QSIM_BLAS_COMPILE_DEFINITIONS)
    target_compile_definitions(blas_compiler_flags INTERFACE ${QSIM_BLAS_COMPILE_DEFINITIONS})
endif()

if(QSIM_BLAS_INCLUDE_DIRS)
    target_include_directories(blas_compiler_flags INTERFACE ${QSIM_BLAS_INCLUDE_DIRS})
endif()

target_link_libraries(blas_compiler_flags INTERFACE ${QSIM_BLAS_LIBRARIES})

# Interface library to specify usage of OpenMP routines
add_library(omp_compiler_flags INTERFACE)

if(ENABLE_OPENMP)
    message(STATUS "OpenMP parallelization enabled")
    target_compile_options(omp_compiler_flags
            INTERFACE
                $<$<BOOL:${APPLE}>:-Xpreprocessor>
                -fopenmp
    )
    target_include_directories(omp_compiler_flags
            INTERFACE $<$<BOOL:${APPLE}>:/opt/homebrew/opt/libomp/include>
    )
    target_link_libraries(omp_compiler_flags
            INTERFACE $<$<BOOL:${APPLE}>:-lomp>
    )
    target_link_directories(omp_compiler_flags
            INTERFACE $<$<BOOL:${APPLE}>:/opt/homebrew/opt/libomp/lib>
    )
    target_link_options(omp_compiler_flags
            INTERFACE $<$<AND:$<BOOL:${UNIX}>,$<NOT:$<BOOL:${APPLE}>>>:-fopenmp>
    )
else()
    message(STATUS "OpenMP parallelization disabled")
endif()
