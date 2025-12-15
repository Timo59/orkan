# CompilerFlags.cmake - Compiler and linking configuration

# Set name for interface library holding compiler flags for the project
set(PROJECT_COMPILER_FLAGS "${PROJECT_NAME}_compiler_flags")

# Interface library for project compiler flags
add_library(${PROJECT_COMPILER_FLAGS} INTERFACE)

target_compile_definitions(${PROJECT_COMPILER_FLAGS}
        INTERFACE
            $<$<BOOL:${APPLE}>:ACCELERATE_NEW_LAPACK>   # Required to use cblas_new
            $<$<BOOL:${APPLE}>:ACCELERATE_LAPACK_ILP64> # __LAPACK_int is 64-bit
            # WARNING: Exposes ILP64 prototypes in header files. OpenBLAS chooses 32-bit vs 64-bit at build time; Ensure
            # OpenBLAS library was built with --INTERFACE64=1
            $<$<AND:$<BOOL:${UNIX}>,$<NOT:$<BOOL:${APPLE}>>>:OPENBLAS_USE64BITINT>  # blasint is 64-bit
)
target_compile_options(${PROJECT_COMPILER_FLAGS}
        INTERFACE
            -O2
            -g
            -fno-strict-aliasing
)
target_link_libraries(${PROJECT_COMPILER_FLAGS}
        INTERFACE $<$<AND:$<BOOL:${UNIX}>,$<NOT:$<BOOL:${APPLE}>>>:-lm>
)

# Interface library to specify usage of BLAS routines
add_library(blas_compiler_flags INTERFACE)

# Use BLAS configuration from Dependencies.cmake
if(QSIM_BLAS_INCLUDE_DIRS)
    target_include_directories(blas_compiler_flags INTERFACE ${QSIM_BLAS_INCLUDE_DIRS})
endif()

target_link_libraries(blas_compiler_flags INTERFACE ${QSIM_BLAS_LIBRARIES})

# Interface library to specify usage of OpenMP routines
add_library(omp_compiler_flags INTERFACE)

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
