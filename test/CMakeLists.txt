add_library(Test STATIC utils.c gateMat.c pauliMat.c finDiff.c)

add_dependencies(Test Math)
target_include_directories(Test
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
)

target_link_libraries(Test PRIVATE Math)