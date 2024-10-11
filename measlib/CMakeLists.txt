set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Meas OBJECT ${SOURCE_DIR}/exact.c)

target_include_directories(Meas
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        PUBLIC  ${INCLUDE_DIR}
)
list(APPEND MEAS_PRIVATE_LINK_LIBRARIES                         # Lists holding all linker libraries needed to compile
        Circ                                                    # Meas
        Gate
        State
)
list(APPEND MEAS_PUBLIC_LINK_LIBRARIES Q_compiler_flags)
list(APPEND install_includes                                    # List of header files file to be included in install
        ${INCLUDE_DIR}/exact.h
)

# Introduce compile time option of API/specification improved gate operations
set(MEASMETH BLAS_OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE MEASMETH PROPERTY STRINGS BLAS, OMP, BLAS_OMP, TEST, OFF)

# Depending on compile time option, build a static library of the improved gate operations and feed it into the declared
# Gatelib
if(MEASMETH STREQUAL BLAS)
elseif(MEASMETH STREQUAL OMP)
elseif(MEASMETH STREQUAL BLAS_OMP)
    target_compile_definitions(Meas PRIVATE EXACT_BLAS_OMP)

    add_library(ExactBlMp STATIC ${SOURCE_DIR}/exactBlMp.c)     # Static library for combined BLAS and OpenMP improved
    target_include_directories(ExactBlMp                        # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
            ${INCLUDE_DIR}
    )
    target_link_libraries(ExactBlMp
            PRIVATE BLAS_compiler_flags
                    Circ
                    Gate
                    OMP_compiler_flags
                    Q_compiler_flags
                    State
    )
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            OMP_compiler_flags
            ExactBlMp
    )
    list(APPEND install_includes ${INCLUDE_DIR}/exactBlMp.h)

elseif (MEASMETH STREQUAL TEST)
else()
endif()
target_link_libraries(Meas
        PRIVATE MEAS_PRIVATE_LINK_LIBRARIES
        PUBLIC  MEAS_PUBLIC_LINK_LIBRARIES
)

install(FILES "exact.h" DESTINATION include)