# CompilerFlags.cmake - Compiler and linking configuration

# Set name for interface library holding compiler flags for the project
set(PROJECT_COMPILER_FLAGS "${PROJECT_NAME}_compiler_flags")

# Interface library for project compiler flags
# Optimization and debug flags:
#   Debug:   -O0 (explicit) + -g (from CMAKE_C/CXX_FLAGS_DEBUG cmake default)
#   Release: -O3 -DNDEBUG   (from CMAKE_C/CXX_FLAGS_RELEASE cmake default) + -march=native
add_library(${PROJECT_COMPILER_FLAGS} INTERFACE)

target_compile_options(${PROJECT_COMPILER_FLAGS}
        INTERFACE
            -fno-strict-aliasing
            $<$<CONFIG:Debug>:-O0>
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
    if(APPLE)
        # Apple Clang does not ship with OpenMP; locate Homebrew's libomp
        # dynamically so the path works on both Apple Silicon and Intel Macs.
        find_program(BREW brew)
        if(BREW)
            execute_process(
                COMMAND ${BREW} --prefix libomp
                OUTPUT_VARIABLE LIBOMP_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            set(OpenMP_C_FLAGS   "-Xpreprocessor -fopenmp")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_C_LIB_NAMES   "omp")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
            list(APPEND CMAKE_PREFIX_PATH "${LIBOMP_PREFIX}")
        else()
            message(FATAL_ERROR "ENABLE_OPENMP=ON but brew not found. Install libomp via Homebrew.")
        endif()
    endif()

    find_package(OpenMP REQUIRED COMPONENTS C CXX)
    message(STATUS "OpenMP parallelization enabled")
    target_link_libraries(omp_compiler_flags INTERFACE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)

    # FindOpenMP does not always propagate the include directory on macOS;
    # add it explicitly from the resolved Homebrew prefix.
    if(APPLE AND LIBOMP_PREFIX)
        target_include_directories(omp_compiler_flags INTERFACE "${LIBOMP_PREFIX}/include")
    endif()
else()
    message(STATUS "OpenMP parallelization disabled")
endif()
