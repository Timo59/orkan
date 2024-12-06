set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Introduce compile time option of API/specification improved gate operations
set(MEASMETH BLAS_OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE MEASMETH PROPERTY STRINGS BLAS OMP BLAS_OMP TEST OFF)

add_library(Meas OBJECT                                                 # Build the object library contributing to qlib
        ${SOURCE_DIR}/exact.c                                           # Source files with meas not yet optimized and
        $<$<STREQUAL:${MEASMETH},BLAS>:${SOURCE_DIR}/exactBlas.c>       # If BLAS, operations optimized with BLAS
        $<$<STREQUAL:${MEASMETH},BLAS_OMP>:${SOURCE_DIR}/exactBlMp.c>   # If BLAS_OMP, meas optimized with BLAS & OMP
        $<$<STREQUAL:${MEASMETH},OLD>:${SOURCE_DIR}/exactOld.c>         # If OLD, unoptimized measurements
        $<$<STREQUAL:${MEASMETH},OMP>:${SOURCE_DIR}/exactOmp.c>         # If OMP, measurements optimized with OpenMP
        $<$<STREQUAL:${MEASMETH},PROFILE>:${SOURCE_DIR}/profile.c>      # If PROFILE, source file, merging all meas
)
target_include_directories(Meas                                         # Include directories for the object library
        PRIVATE ${CMAKE_SOURCE_DIR}/include                             # qTypes only populates this include directories
        PUBLIC  ${CMAKE_CURRENT_BINARY_DIR}                             # Include files from measlib are propagated
)
target_link_libraries(Meas                                              # Link libraries depending on the Cache string
        PRIVATE State
                Gate
                Circ
                $<$<IN_LIST:${MEASMETH},OLD;OMP;PROFILE>:Math>
        PUBLIC  Q_compiler_flags
                $<$<IN_LIST:${MEASMETH},BLAS;BLAS_OMP;PROFILE>:BLAS_compiler_flags>
                $<$<IN_LIST:${MEASMETH},OMP;BLAS_OMP;PROFILE>:OMP_compiler_flags>
)

if (${MEASMETH} STREQUAL BLAS)
    set(INCLUDE_VECLIB "#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${MEASMETH} STREQUAL BLAS_OMP)
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${MEASMETH} STREQUAL OMP)
    set(INCLUDE_LINALG "#ifndef LINALG_H\n#include \"linalg.h\"\n#endif\n")
    set(INCLUDE_OMP "\n#ifndef __OMP_H\n#include <omp.h>\n#endif")
elseif (${MEASMETH} STREQUAL PROFILE)
    set(INCLUDE_LINALG "#ifndef LINALG_H\n#include \"linalg.h\"\n#endif\n")
    set(INCLUDE_OMP "\n#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif\n")
    set(PAULI_EXPVAL_PROFILE
            "\ncplx_t expValPauli_blas(const state_t* state, const pauliOp_t* op, ptrStr applyStr);
             \ncplx_t expValPauli_blas_omp(const state_t* state, const pauliOp_t* op, ptrStr applyStr);
             \ncplx_t expValPauli_old(const state_t* state, const pauliOp_t* op, ptrStr applyStr);
             \ncplx_t expValPauli_omp(const state_t* state, const pauliOp_t* op, ptrStr applyStr);\n"
    )
    set(PAULI_EXPVALOBS_PROFILE
            "\ndouble expValObsPauli_blas(const state_t* state, const pauliObs_t* obs, ptrStr applyStr);
             \ndouble expValObsPauli_blas_omp(const state_t* state, const pauliObs_t* obs, ptrStr applyStr);
             \ndouble expValObsPauli_old(const state_t* state, const pauliObs_t* obs, ptrStr applyStr);
             \ndouble expValObsPauli_omp(const state_t* state, const pauliObs_t* obs, ptrStr applyStr);\n"
    )
    set(PAULI_GRAD_PROFILE
            "\ndouble* gradPQC_blas(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                     depth_t circdepth);
             \ndouble* gradPQC_blas_omp(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                         depth_t circdepth);
             \ndouble* gradPQC_omp(const state_t* state, const double par[], const obs_t* obs, const obs_t *evoOps[],
                    depth_t circdepth);\n"
    )
else ()
    set(INCLUDE_LINALG "#ifndef LINALG_H\n#include \"linalg.h\"\n#endif")
endif ()

configure_file(${INCLUDE_DIR}/exact.h.in exact.h @ONLY)

install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/exact.h
        DESTINATION include/qlib)