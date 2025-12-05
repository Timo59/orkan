# Find unity.h
find_path(UNITY_DIR
        NAMES unity.h
        PATHS $ENV{HOME}/Code/Unity/src
)
if (NOT UNITY_DIR)
    message(FATAL_ERROR "Unity source directory not found.")
endif ()

# Source files for test_q
set(TEST_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/test_state.c"
        "${UNITY_DIR}/unity.c"
)

add_executable(test_q ${TEST_SOURCES})
target_compile_definitions(test_q
        PRIVATE
            UNITY_INCLUDE_DOUBLE
            UNITY_SUPPORT_64
)
target_include_directories(test_q
        PRIVATE
            ${UNITY_DIR}
            "${CMAKE_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(test_q PRIVATE q)