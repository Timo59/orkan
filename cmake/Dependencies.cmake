# Dependencies.cmake - Handle all external dependencies

# Option to use system OpenBLAS (LP64; only used by tests/benchmarks for small matrices)
option(USE_SYSTEM_OPENBLAS "Use system OpenBLAS instead of bundled" OFF)

# BLAS/LAPACK configuration
if(APPLE)
    # Apple Accelerate framework has guaranteed ILP64 support with ACCELERATE_LAPACK_ILP64
    message(STATUS "Using Apple Accelerate framework (ILP64 native)")
    set(QSIM_BLAS_COMPILE_DEFINITIONS
            ACCELERATE_NEW_LAPACK # Required to use cblas_new
            ACCELERATE_LAPACK_ILP64  # __LAPACK_int is 64-bit
    )
    set(QSIM_BLAS_LIBRARIES "-framework Accelerate")
    set(QSIM_BLAS_INCLUDE_DIRS "")

elseif(UNIX)
    # Linux: Bundle OpenBLAS with ILP64 by default
    if(USE_SYSTEM_OPENBLAS)
        message(STATUS "Using system OpenBLAS (LP64)")

        # Try to find system OpenBLAS
        find_library(OPENBLAS_LIB
            NAMES openblas
            PATHS /usr/lib /usr/local/lib /opt/OpenBLAS/lib
            DOC "OpenBLAS library"
        )

        if(NOT OPENBLAS_LIB)
            message(FATAL_ERROR
                "System OpenBLAS not found. Either:\n"
                "  1. Install: sudo apt install libopenblas-dev\n"
                "  2. Build with: cmake -DUSE_SYSTEM_OPENBLAS=OFF .. (uses bundled ILP64 OpenBLAS)")
        endif()

        message(STATUS "Found OpenBLAS: ${OPENBLAS_LIB}")
        set(QSIM_BLAS_COMPILE_DEFINITIONS "")
        set(QSIM_BLAS_LIBRARIES ${OPENBLAS_LIB})
        set(QSIM_BLAS_INCLUDE_DIRS "/usr/include" "/usr/local/include")

    else()
        # Use bundled OpenBLAS with ILP64 support
        message(STATUS "Using bundled OpenBLAS with ILP64 support")
        message(STATUS "Note: First build will take a few minutes to compile OpenBLAS")

        include(FetchContent)
        FetchContent_Declare(
            openblas
            URL https://github.com/xianyi/OpenBLAS/releases/download/v0.3.27/OpenBLAS-0.3.27.tar.gz
            URL_HASH SHA256=aa2d68b1564fe2b13bc292672608e9cdeeeb6dc34995512e65c3b10f4599e897
        )

        # OpenBLAS build options
        set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
        set(BUILD_STATIC_LIBS OFF CACHE BOOL "Build static libraries" FORCE)
        set(INTERFACE64 1 CACHE STRING "Use 64-bit integers (ILP64)" FORCE)
        set(DYNAMIC_ARCH OFF CACHE BOOL "Build for multiple CPU architectures" FORCE)
        set(BUILD_TESTING OFF CACHE BOOL "Disable OpenBLAS tests" FORCE)
        set(BUILD_WITHOUT_LAPACK OFF CACHE BOOL "Include LAPACK" FORCE)

        FetchContent_MakeAvailable(openblas)

        # OpenBLAS creates the target 'openblas'
        set(QSIM_BLAS_COMPILE_DEFINITIONS OPENBLAS_USE64BITINT)
        set(QSIM_BLAS_LIBRARIES openblas)
        # Get the include directory from the fetched content
        FetchContent_GetProperties(openblas SOURCE_DIR OPENBLAS_SOURCE_DIR)
        set(QSIM_BLAS_INCLUDE_DIRS "${OPENBLAS_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/generated")
    endif()

else()
    message(FATAL_ERROR "Unsupported platform: ${CMAKE_SYSTEM_NAME}")
endif()

message(STATUS "BLAS configuration complete")
