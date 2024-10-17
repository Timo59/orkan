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
target_compile_definitions(Circ
        PUBLIC  $<$<STREQUAL:${CIRCMETH},BLAS>:PAULI_BLAS>              # Add the compile definition depending on Cache
                $<$<STREQUAL:${CIRCMETH},BLAS_OMP>:PAULI_BLAS_OMP>      # string
                $<$<STREQUAL:${CIRCMETH},OMP>:PAULI_OMP>
                $<$<STREQUAL:${CIRCMETH},PROFILE>:PAULI_PROFILE>
)
target_include_directories(Circ                                         # Include directories for the object library
        PRIVATE "${CMAKE_SOURCE_DIR}/include"                           # qTypes only populates this include directories
        PUBLIC  ${INCLUDE_DIR}                                          # Include files from circlib are propagated
)
target_link_libraries(Circ                                              # Link libraries depending on the Cache string
        PRIVATE State
                Gate
        PUBLIC  BLAS_compiler_flags                     # needs to be included as long as lcupqg is part of circlib.c
                Q_compiler_flags
                $<$<IN_LIST:${CIRCMETH},OLD;OMP;PROFILE>:Math>
                $<$<IN_LIST:${CIRCMETH},BLAS_OMP;OMP;PROFILE>:OMP_compiler_flags>
)

install(FILES       ${INCLUDE_DIR}/circlib.h
                    ${INCLUDE_DIR}/diag.h
                    ${INCLUDE_DIR}/pauli.h
        DESTINATION include
)
install(FILES $<$<STREQUAL:${CIRCMETH},PROFILE>:${INCLUDE_DIR}/profile.h>   # Install file depending on the Cache string
        DESTINATION include
        RENAME      profile_pauli.h
)