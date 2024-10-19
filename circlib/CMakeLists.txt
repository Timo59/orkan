set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Introduce compile time option of API/specification improved gate operations
set(CIRCMETH BLAS_OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE CIRCMETH PROPERTY STRINGS BLAS OMP BLAS_OMP TEST OFF)

add_library(Circ OBJECT                                                 # Build the object library contributing to qlib
        ${SOURCE_DIR}/circlib.c                                         # Source files with operations not yet optimized
        ${SOURCE_DIR}/diag.c                                            # and
        ${SOURCE_DIR}/pauli.c
        $<$<STREQUAL:${CIRCMETH},BLAS>:${SOURCE_DIR}/pauliBlas.c>       # If BLAS, operations optimized with BLAS
        $<$<STREQUAL:${CIRCMETH},BLAS_OMP>:${SOURCE_DIR}/pauliBlMp.c>   # If BLAS_OMP, ops optimized with BLAS & OMP
        $<$<STREQUAL:${CIRCMETH},OLD>:${SOURCE_DIR}/pauliOld.c>         # If OLD, unoptimized operations
        $<$<STREQUAL:${CIRCMETH},OMP>:${SOURCE_DIR}/pauliOmp.c>         # If OMP, operations optimized with OpenMP
        $<$<STREQUAL:${CIRCMETH},PROFILE>:${SOURCE_DIR}/profile.c>      # If PROFILE, source file, merging all ops
)
target_include_directories(Circ                                         # Include directories for the object library
        PRIVATE ${CMAKE_SOURCE_DIR}/include                             # qTypes only populates this include directories
        PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include                     # Include files from circlib are propagated
                ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(Circ                                              # Link libraries depending on the Cache string
        PRIVATE State
                Gate
        PUBLIC  BLAS_compiler_flags                     # needs to be included as long as lcupqg is part of circlib.c
                Q_compiler_flags
                $<$<IN_LIST:${CIRCMETH},OLD;OMP;PROFILE>:Math>
                $<$<IN_LIST:${CIRCMETH},BLAS_OMP;OMP;PROFILE>:OMP_compiler_flags>
)

if (${CIRCMETH} STREQUAL BLAS)
    set(INCLUDE_LINALG "")
    set(INCLUDE_OMP "")
    set(INCLUDE_VECLIB "#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${CIRCMETH} STREQUAL BLAS_OMP)
    set(INCLUDE_LINALG "")
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")
elseif (${CIRCMETH} STREQUAL OMP)
    set(INCLUDE_LINALG "#ifndef LINALG_H\n#include \"linalg.h\"\n#endif\n")
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif")
    set(INCLUDE_VECLIB "")
elseif (${CIRCMETH} STREQUAL PROFILE)
    set(INCLUDE_LINALG "#ifndef LINALG_H\n#include \"linalg.h\"\n#endif\n")
    set(INCLUDE_OMP "#ifndef __OMP_H\n#include <omp.h>\n#endif\n")
    set(INCLUDE_VECLIB "\n#ifndef __VECLIB__\n#include <vecLib/vecLib.h>\n#endif")

    set(PAULI_OP_PROFILE "\nvoid applyOpPauli_blas(state_t* state, const pauliOp_t* op);
                          \nvoid applyOpPauli_blas_omp(state_t* state, const pauliOp_t* op);
                          \nvoid applyOpPauli_omp(state_t* state, const pauliOp_t* op);")
    set(PAULI_OBS_PROFILE "\nvoid applyObsPauli_blas(state_t* state, const pauliObs_t* obs);
                           \nvoid applyObsPauli_blas_omp(state_t* state, const pauliObs_t* obs);
                           \nvoid applyObsPauli_omp(state_t* state, const pauliObs_t* obs);")
    set(PAULI_EVO_PROFILE "\nvoid evolvePauliStr_blas(state_t* state, const pauli_t paulistr[], double angle);
                           \nvoid evolvePauliStr_blas_omp(state_t* state, const pauli_t paulistr[], double angle);
                           \nvoid evolvePauliStr_omp(state_t* state, const pauli_t paulistr[], double angle);")
endif ()

configure_file(${INCLUDE_DIR}/pauli.h.in pauli.h @ONLY)
install(FILES       ${INCLUDE_DIR}/circlib.h
                    ${INCLUDE_DIR}/diag.h
                    ${CMAKE_CURRENT_BINARY_DIR}/pauli.h
        DESTINATION include
)
