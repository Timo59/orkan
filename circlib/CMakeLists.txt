set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Circ OBJECT
        ${SOURCE_DIR}/diag.c
        ${SOURCE_DIR}/pauli.c
        ${SOURCE_DIR}/circlib.c
)
target_include_directories(Circ
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        PUBLIC  ${INCLUDE_DIR}
)
list(APPEND CIRC_PRIVATE_LINK_LIBRARIES                         # Lists holding all linker libraries needed to compile
        Gate                                                    # Circ
        State
)
list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
        BLAS_compiler_flags                                     # needs to be included as long as lcupqg is part of
                                                                # of circlib.c
        Q_compiler_flags
)
list(APPEND install_includes                                    # List of header files file to be included in install
        ${INCLUDE_DIR}/circlib.h
        ${INCLUDE_DIR}/diag.h
        ${INCLUDE_DIR}/pauli.h
)

# Introduce compile time option of API/specification improved gate operations
set(CIRCMETH OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE CIRCMETH PROPERTY STRINGS BLAS, OMP, BLAS_OMP, TEST, OFF)

# Depending on compile time option, build a static library of the improved gate operations and feed it into the declared
# Gatelib
if(CIRCMETH STREQUAL BLAS)
    target_compile_definitions(Circ PRIVATE PAULI_BLAS)

    add_library(PauliBlas STATIC ${SOURCE_DIR}/pauliBlas.c)     # Static library for combined BLAS and OpenMP improved
    target_include_directories(PauliBlas                        # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
            ${INCLUDE_DIR}
    )
    target_link_libraries(PauliBlas
            PRIVATE BLAS_compiler_flags
                    Gate
                    Q_compiler_flags
                    State
    )
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            PauliBlas
    )
    list(APPEND install_includes ${INCLUDE_DIR}/pauliBlas.h)

elseif(CIRCMETH STREQUAL OMP)
    target_compile_definitions(Circ PRIVATE PAULI_OMP)

    add_library(PauliOmp STATIC ${SOURCE_DIR}/pauliOmp.c)       # Static library for combined BLAS and OpenMP improved
    target_include_directories(PauliOmp                         # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
            ${INCLUDE_DIR}
    )
    target_link_libraries(PauliOmp
            PRIVATE Gate
                Math
                OMP_compiler_flags
                Q_compiler_flags
                State
    )
    list(APPEND CIRC_PRIVATE_LINK_LIBRARIES Math)
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
            OMP_compiler_flags
            PauliOmp
    )
    list(APPEND ${INCLUDE_DIR}/pauliOmp.h)

elseif(CIRCMETH STREQUAL BLAS_OMP)
    target_compile_definitions(Circ PRIVATE PAULI_BLAS_OMP)

    add_library(PauliBlMp STATIC ${SOURCE_DIR}/pauliBlMp.c)     # Static library for combined BLAS and OpenMP improved
    target_include_directories(PauliBlMp                        # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(PauliBlMp
            PRIVATE BLAS_compiler_flags
                    Gate
                    OMP_compiler_flags
                    Q_compiler_flags
                    State
    )
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            OMP_compiler_flags
            PauliBlMp
    )
    list(APPEND install_includes ${INCLUDE_DIR}/pauliBlMp.h)
elseif(CIRCMETH STREQUAL TEST)
    target_compile_definitions(Circ PRIVATE TEST)

    add_library(PauliTest STATIC                                # Static library for test inherits functions from all
            ${SOURCE_DIR}/pauliBlas.c                           # improved gate operations' source files
            ${SOURCE_DIR}/pauliOmp.c
            ${SOURCE_DIR}/pauliBlMp.c
            ${SOURCE_DIR}/pauliOld.c
    )
    target_include_directories(PauliTest
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
            ${INCLUDE_DIR}
    )
    target_link_libraries(PauliTest
            PRIVATE BLAS_compiler_flags
            Gate
            Math
            OMP_compiler_flags
            Q_compiler_flags
            State
    )

    list(APPEND CIRC_PRIVATE_LINK_LIBRARIES Math)
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            OMP_compiler_flags
            PauliTest
    )
    list(APPEND install_includes
            ${INCLUDE_DIR}/pauliBlas.h
            ${INCLUDE_DIR}/pauliBlMp.h
            ${INCLUDE_DIR}/pauliOld.h
            ${INCLUDE_DIR}/pauliOmp.h
    )

else()
    add_library(PauliOld STATIC ${SOURCE_DIR}/pauliOld.c)       # Static library for combined BLAS and OpenMP improved
    target_include_directories(PauliOld                         # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
            ${INCLUDE_DIR}
    )
    target_link_libraries(PauliOld
            PRIVATE Gate
                    Math
                    Q_compiler_flags
                    State
    )
    list(APPEND CIRC_PRIVATE_LINK_LIBRARIES Math)
    list(APPEND CIRC_PUBLIC_LINK_LIBRARIES PauliOld)
    list(APPEND install_includes ${INCLUDE_DIR}/pauliOld.h)
endif()

target_link_libraries(Circ
        PRIVATE ${CIRC_PRIVATE_LINK_LIBRARIES}
        PUBLIC  ${CIRC_PUBLIC_LINK_LIBRARIES}
)

install(FILES "circlib.h" "diag.h" "pauli.h" DESTINATION include)