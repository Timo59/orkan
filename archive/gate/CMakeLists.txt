set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Gate OBJECT                                                 # Build the object library contributing to qlib
        ${SOURCE_DIR}/gatelib.c                                         # Source file with gates not yet optimized and
        $<$<STREQUAL:${GATEMETH},BLAS>:${SOURCE_DIR}/qhipBLAS.c>        # If BLAS, gates optimized with BLAS
        $<$<STREQUAL:${GATEMETH},BLAS_OMP>:${SOURCE_DIR}/qhipBlMp.c>    # If BLAS_OMP, gates optimized with BLAS & OMP
        $<$<STREQUAL:${GATEMETH},GCD>:${SOURCE_DIR}/qhipGcd.c>          # If GCD, gates optimized with GCD
        $<$<STREQUAL:${GATEMETH},OLD>:${SOURCE_DIR}/qhipster.c>         # If OLD, unoptimized gates
        $<$<STREQUAL:${GATEMETH},OMP>:${SOURCE_DIR}/qhipOmp.c>          # If OMP, gates optimized with OpenMP
        $<$<STREQUAL:${GATEMETH},PROFILE>:${SOURCE_DIR}/profile.c>      # If PROFILE, src file, merging all gates
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

if (${GATEMETH} STREQUAL PROFILE)
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

configure_file(${INCLUDE_DIR}/qhipster.h.in qhipster.h @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/qhipster.h DESTINATION include/qlib)