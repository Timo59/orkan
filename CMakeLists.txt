cmake_minimum_required(VERSION 3.27)

# Set the name of the project and store it to CMAKE_PROJECT_NAME. Set PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR
project(QSim
        VERSION     1.0.0   # Sets PROJECT_VERSION and <PROJECT_NAME>_VERSION
        DESCRIPTION "Quantum Computer simulator" # Sets PROJECT_DESCRIPTION and <PROJECT_NAME>_DESCRIPTION
        LANGUAGES   C CXX   # Programming languages to build the project
)

set(CMAKE_C_STANDARD 17)

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Options
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" ON)
option(ENABLE_ASAN "Enable AddressSanitizer for memory debugging" OFF)

if(ENABLE_ASAN)
    message(STATUS "AddressSanitizer enabled")
    set(SANITIZE_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SANITIZE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZE_FLAGS}")
    # Optional: add for linker flags too (needed on some Linux/macOS setups)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZE_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SANITIZE_FLAGS}")
endif()

# Platform detection and compiler setup
include(cmake/PlatformConfig.cmake)

# Find dependencies
include(cmake/Dependencies.cmake)

# Create compiler flags interface library
include(cmake/CompilerFlags.cmake)

# Build main library
add_subdirectory(src)

if (BUILD_TESTS AND EXISTS "${CMAKE_SOURCE_DIR}/test")
    enable_testing()
    add_subdirectory(test)

    # Require tests to pass before install
    add_custom_target(check
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    )

    add_dependencies(install check)
endif ()

if(BUILD_EXAMPLES AND EXISTS "${CMAKE_SOURCE_DIR}/examples")
    add_subdirectory(examples)
endif()

#add_subdirectory(internal)
#target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)                       # include double has to be called
#                                                                                    # manually when building unity lib
#add_library(BLAS_compiler_flags INTERFACE)                                          # Interface library specifying usage
#target_include_directories(BLAS_compiler_flags                                      # of BLAS routines
#     INTERFACE   $<${LINUX}:/usr/local/include>
#)
#target_link_directories(BLAS_compiler_flags
#     INTERFACE   $<${LINUX}:/usr/local/lib>
#)
#target_link_libraries(BLAS_compiler_flags INTERFACE
#        $<${MACOS}:-lblas>
#        $<${LINUX}:-lopenblas>
#        # $<${LINUX}:-lpthread>
#        $<${MACOS}:-llapack>
#)
#
#add_library(OMP_compiler_flags INTERFACE)                                           # Interface library specifying usage
#target_compile_options(OMP_compiler_flags INTERFACE                                 # requirements for dependents of
#        $<${MACOS}:-Xpreprocessor>                                                  # the OpenMP API
#        -fopenmp
#)
#target_include_directories(OMP_compiler_flags INTERFACE $<${MACOS}:/opt/homebrew/opt/libomp/include>)
#target_link_libraries(OMP_compiler_flags INTERFACE $<${MACOS}:-lomp>)
#target_link_options(OMP_compiler_flags INTERFACE $<${LINUX}:-fopenmp>)
#target_link_directories(OMP_compiler_flags INTERFACE $<${MACOS}:/opt/homebrew/opt/libomp/lib>)
#
#add_library(Q_compiler_flags INTERFACE)                                             # Interface library specifying usage
#target_compile_options(Q_compiler_flags INTERFACE -O2 -g)                           # requirements for all Q dependents
#
#add_library(Test_flags INTERFACE)                                                   # Interface library specifying usage
#target_include_directories(Test_flags                                               # requirements for unit tests
#        INTERFACE   ${CMAKE_SOURCE_DIR}/include
#                    ${CMAKE_SOURCE_DIR}/internal/include
#)
#target_link_libraries(Test_flags INTERFACE
#        Internal
#        unity
#)
#
#add_subdirectory(state)
#add_subdirectory(gate)
#add_subdirectory(circ)
#add_subdirectory(meas)
#
#add_test(NAME testGate COMMAND gateTest)
#add_test(NAME testCirc COMMAND circTest)
#add_test(NAME testMeas COMMAND measTest)
#
#add_library(${PROJECT_NAME} SHARED
#        $<TARGET_OBJECTS:State>
#        $<TARGET_OBJECTS:Gate>
#        $<TARGET_OBJECTS:Circ>
#        $<TARGET_OBJECTS:Meas>
#)
#target_link_libraries(${PROJECT_NAME}
#        PUBLIC  State
#                Gate
#                Circ
#                Meas
#)
#
#install(TARGETS Q
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib
#        RUNTIME DESTINATION bin
#        INCLUDES DESTINATION include/qlib
#)
#install(FILES       include/qLib.h
#                    include/qTypes.h
#        DESTINATION include/qlib)
#
#
#add_dependencies(Q
#        gateTest
#        circTest
#        measTest
#)
#add_custom_command(
#        TARGET Q
#        COMMENT "Run tests on qlib"
#        PRE_BUILD
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
#)
