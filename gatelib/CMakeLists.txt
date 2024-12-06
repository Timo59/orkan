set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Introduce compile time option of API/specification improved gate operations
set(GATEMETH OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE GATEMETH PROPERTY STRINGS BLAS BLAS_OMP GCD OMP PROFILE OLD)

add_library(Gate OBJECT                                                 # Build the object library contributing to qlib
        ${SOURCE_DIR}/gatelib.c                                         # Source file with gates not yet optimized and
        $<$<STREQUAL:${GATEMETH},BLAS>:${SOURCE_DIR}/qhipBLAS.c>        # If BLAS, gates optimized with BLAS
        $<$<STREQUAL:${GATEMETH},BLAS_OMP>:${SOURCE_DIR}/qhipBlMp.c>    # If BLAS_OMP, gates optimized with BLAS & OMP
        $<$<STREQUAL:${GATEMETH},GCD>:${SOURCE_DIR}/qhipGcd.c>          # If GCD, gates optimized with GCD
        $<$<STREQUAL:${GATEMETH},OLD>:${SOURCE_DIR}/qhipster.c>         # If OLD, unoptimized gates
        $<$<STREQUAL:${GATEMETH},OMP>:${SOURCE_DIR}/qhipOmp.c>          # If OMP, gates optimized with OpenMP
        $<$<STREQUAL:${GATEMETH},PROFILE>:${SOURCE_DIR}/profile.c>      # If PROFILE, source file, merging all gates
)
target_include_directories(Gate                                 # Include directories for the object library
        PRIVATE "${CMAKE_SOURCE_DIR}/include"                   # qTypes only populates the include directories for Gate
        PUBLIC  ${CMAKE_CURRENT_BINARY_DIR}                     # Include files from gatelib are propagated
)
target_link_libraries(Gate                                      # Link libraries depending on the Cache string
        PRIVATE State
        PUBLIC  Q_compiler_flags
                $<$<IN_LIST:${GATEMETH},BLAS;BLAS_OMP;PROFILE>:BLAS_compiler_flags>
                $<$<IN_LIST:${GATEMETH},BLAS_OMP;OMP;PROFILE>:OMP_compiler_flags>
)

if (${GATEMETH} STREQUAL BLAS)
    set(INCLUDE_GCD "")
    set(INCLUDE_OMP "")
    set(INCLUDE_VECLIB "#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${GATEMETH} STREQUAL BLAS_OMP)
    set(INCLUDE_GCD "")
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${GATEMETH} STREQUAL GCD)
    set(INCLUDE_GCD "#ifndef __DISPATCH_PUBLIC__\n#include <dispatch/dispatch.h>\n#endif")
    set(INCLUDE_OMP "")
    set(INCLUDE_VECLIB "")
elseif (${GATEMETH} STREQUAL OMP)
    set(INCLUDE_GCD "")
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif")
    set(INCLUDE_VECLIB "")
elseif (${GATEMETH} STREQUAL PROFILE)
    set(INCLUDE_GCD "#ifndef __DISPATCH_PUBLIC__\n#include <dispatch/dispatch.h>\n#endif\n")
    set(INCLUDE_OMP "\n#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")

    set(PAULI_X_PROFILE "\nvoid applyX_blas(state_t* state, qubit_t qubit);
                         \nvoid applyX_blas_omp(state_t* state, qubit_t qubit);
                         \nvoid applyX_gcd(state_t* state, qubit_t qubit);
                         \nvoid applyX_omp(state_t* state, qubit_t qubit);\n")
    set(PAULI_Y_PROFILE "\nvoid applyY_blas(state_t* state, qubit_t qubit);
                         \nvoid applyY_blas_omp(state_t* state, qubit_t qubit);
                         \nvoid applyY_gcd(state_t* state, qubit_t qubit);
                         \nvoid applyY_omp(state_t* state, qubit_t qubit);\n")
    set(PAULI_Z_PROFILE "\nvoid applyZ_blas(state_t* state, qubit_t qubit);
                         \nvoid applyZ_blas_omp(state_t* state, qubit_t qubit);
                         \nvoid applyZ_gcd(state_t* state, qubit_t qubit);
                         \nvoid applyZ_omp(state_t* state, qubit_t qubit);\n")
endif ()

configure_file(${INCLUDE_DIR}/gatelib.h.in gatelib.h @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gatelib.h DESTINATION include/qlib)