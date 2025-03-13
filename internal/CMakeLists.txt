set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(Internal OBJECT
        ${SOURCE_DIR}/finDiff.c
        ${SOURCE_DIR}/gateMat.c
        ${SOURCE_DIR}/linalg.c
        ${SOURCE_DIR}/pauliMat.c
        ${SOURCE_DIR}/utils.c
)
target_include_directories(Internal
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PUBLIC  ${INCLUDE_DIR}
)