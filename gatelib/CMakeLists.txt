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
target_compile_definitions(Gate                                         # Add the compile definition depending on Cache
        PUBLIC  $<$<STREQUAL:${GATEMETH},BLAS>:GATE_BLAS>               # string
                $<$<STREQUAL:${GATEMETH},BLAS_OMP>:GATE_BLAS_OMP>
                $<$<STREQUAL:${GATEMETH},GCD>:GATE_GCD>
                $<$<STREQUAL:${GATEMETH},OMP>:GATE_OMP>
                $<$<STREQUAL:${GATEMETH},PROFILE>:GATE_PROFILE>
)
target_include_directories(Gate                                 # Include directories for the object library
        PRIVATE "${CMAKE_SOURCE_DIR}/include"                   # qTypes only populates the include directories for Gate
        PUBLIC  ${INCLUDE_DIR}                                  # Include files from gatelib are propagated
)
target_link_libraries(Gate                                      # Link libraries depending on the Cache string
        PRIVATE State
        PUBLIC  Q_compiler_flags
                $<$<IN_LIST:${GATEMETH},BLAS;BLAS_OMP;PROFILE>:BLAS_compiler_flags>
                $<$<IN_LIST:${GATEMETH},BLAS_OMP;OMP;PROFILE>:OMP_compiler_flags>
)

install(FILES ${INCLUDE_DIR}/gatelib.h DESTINATION include)
install(FILES $<$<STREQUAL:${GATEMETH},PROFILE>:${INCLUDE_DIR}/profile.h>   # Install file depending on the Cache string
        DESTINATION include
        RENAME      profile_gate.h
)
