set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Gate OBJECT                                         # Object library of gate operation objects
        ${SOURCE_DIR}/gatelib.c
        ${SOURCE_DIR}/qhipster.c
)
target_include_directories(Gate                                 # Include directories for the object library
        PRIVATE "${CMAKE_SOURCE_DIR}/include"                   # qTypes only populates the include directories for Gate
        PUBLIC  ${INCLUDE_DIR}                                  # Include files from gatelib are propagated
)
list(APPEND GATE_PUBLIC_LINK_LIBRARIES Q_compiler_flags)
list(APPEND install_includes                                    # Set a list holding include files for install
        gatelib.h
        qhipster.h
)

# Introduce compile time option of API/specification improved gate operations
set(GATEMETH OMP CACHE STRING "Select the method to apply quantum gates")
set_property(CACHE GATEMETH PROPERTY STRINGS BLAS, OMP, GCD, BLAS_OMP, TEST, OFF)

# Depending on compile time option, build a static library of the improved gate operations and feed it into the declared
# Gatelib
if(GATEMETH STREQUAL BLAS)
    target_compile_definitions(Gate PRIVATE GATE_BLAS)

    add_library(GateBlas STATIC ${SOURCE_DIR}/qhipBlas.c)       # Static library for BLAS improved gate operations
    target_include_directories(GateBlas
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(GateBlas
            PRIVATE BLAS_compiler_flags
                    Q_compiler_flags
                    State
    )
    list(APPEND GATE_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            GateBlas
    )
    list(APPEND install_includes qhipBlas.h)

elseif(GATEMETH STREQUAL OMP)
    target_compile_definitions(Gate PRIVATE GATE_OMP)

    add_library(GateOmp STATIC ${SOURCE_DIR}/qhipOmp.c)         # Static library for OpenMP improved gate operations
    target_include_directories(GateOmp
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(GateOmp
            PRIVATE OMP_compiler_flags
                    Q_compiler_flags
                    State
    )

    list(APPEND GATE_PUBLIC_LINK_LIBRARIES
            GateOmp
            OMP_compiler_flags
    )
    list(APPEND install_includes qhipOmp.h)

elseif(GATEMETH STREQUAL GCD)
    target_compile_definitions(Gate PRIVATE GATE_GCD)

    add_library(GateGcd STATIC ${SOURCE_DIR}/qhipGcd.c)         # Static library for Grand Central Dispatch improved
    target_include_directories(GateGcd                          # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(GateGcd PRIVATE State)

    list(APPEND GATE_PUBLIC_LINK_LIBRARIES GateGcd)
    list(APPEND install_includes qhipGcd.h)

elseif(GATEMETH STREQUAL BLAS_OMP)
    target_compile_definitions(Gate PRIVATE GATE_BLAS_OMP)

    add_library(GateBlMp STATIC ${SOURCE_DIR}/qhipBlMp.c)       # Static library for combined BLAS and OpenMP improved
    target_include_directories(GateBlMp                         # gate operations
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(GateBlMp
            PRIVATE BLAS_compiler_flags
                    OMP_compiler_flags
                    Q_compiler_flags
                    State
    )

    list(APPEND GATE_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            GateBlMp
            OMP_compiler_flags
    )
    list(APPEND install_includes qhipBlMP.h)

elseif(GATEMETH STREQUAL TEST)
    target_compile_definitions(Gate PRIVATE TEST)

    add_library(GateTest STATIC                                 # Static library for test inherits functions from all
            ${SOURCE_DIR}/qhipBlas.c                            # improved gate operations' source files
            ${SOURCE_DIR}/qhipOmp.c
            ${SOURCE_DIR}/qhipGcd.c
            ${SOURCE_DIR}/qhipBlMp.c
    )
    target_include_directories(GateTest
            PRIVATE "${CMAKE_SOURCE_DIR}/include"
                    ${INCLUDE_DIR}
    )
    target_link_libraries(GateTest
            PRIVATE BLAS_compiler_flags
                    OMP_compiler_flags
                    Q_compiler_flags
                    State
    )

    list(APPEND GATE_PUBLIC_LINK_LIBRARIES
            BLAS_compiler_flags
            GateTest
            OMP_compiler_flags
    )
    list(APPEND install_includes
            ${INCLUDE_DIR}/qhipBlas.h
            ${INCLUDE_DIR}/qhipOmp.h
            ${INCLUDE_DIR}/qhipGcd.h
            ${INCLUDE_DIR}/qhipBlMp.h
    )
endif()

target_link_libraries(Gate
        PRIVATE State
        PUBLIC ${GATE_PUBLIC_LINK_LIBRARIES}
)

install(FILES ${install_includes} DESTINATION include)