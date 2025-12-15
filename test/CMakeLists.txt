# Initialize submodule if not present
if (NOT EXISTS PATHS "${CMAKE_SOURCE_DIR}/extern/Unity/src/unity.h")
    find_package(Git QUIET)
    if(GIT_FOUND)
        message(STATUS "Initializing Unity submodule...")
        execute_process(
                COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive extern/Unity
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT
        )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}")
        endif()
    else()
        message(FATAL_ERROR "Git not found. Cannot initialize Unity submodule")
    endif()
endif()

# Find unity.h
find_path(UNITY_SRC
        NAMES unity.h
        PATHS "${CMAKE_SOURCE_DIR}/extern/Unity/src"
        NO_DEFAULT_PATH
)
if (NOT UNITY_SRC)
    message(FATAL_ERROR "Unity source directory not found.")
endif ()

# Source files for test_q
set(TEST_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/test.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/test_state.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/test_qhipster.c"
        "${UNITY_SRC}/unity.c"
)

add_executable(test_q ${TEST_SOURCES})
target_compile_definitions(test_q
        PRIVATE
            UNITY_INCLUDE_DOUBLE
            UNITY_SUPPORT_64
)
target_include_directories(test_q
        PRIVATE
            ${UNITY_SRC}
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_q PRIVATE q)