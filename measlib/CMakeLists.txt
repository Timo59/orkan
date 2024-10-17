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
target_compile_definitions(Meas                                         # Add the compile definition depending on Cache
        PUBLIC  $<$<STREQUAL:${MEASMETH},BLAS>:EXACT_BLAS>              # string
                $<$<STREQUAL:${MEASMETH},BLAS_OMP>:EXACT_BLAS_OMP>
                $<$<STREQUAL:${MEASMETH},OMP>:EXACT_OMP>
                $<$<STREQUAL:${MEASMETH},PROFILE>:EXACT_PROFILE>
)
target_include_directories(Meas                                         # Include directories for the object library
        PRIVATE "${CMAKE_SOURCE_DIR}/include"                           # qTypes only populates this include directories
        PUBLIC  ${INCLUDE_DIR}                                          # Include files from measlib are propagated
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

install(FILES       ${INCLUDE_DIR}/exact.h
        DESTINATION include)
install(FILES $<$<STREQUAL:${MEASMETH},PROFILE>:${INCLUDE_DIR}/profile.h> # Install file depending on the Cache string
        DESTINATION include
        RENAME      profile_exact.h
)