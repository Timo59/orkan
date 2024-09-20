add_library(Test STATIC
        utils.c
        $<TARGET_OBJECTS:Math>
        gateMat.c
        pauliMat.c
        finDiff.c)

add_dependencies(Test Math)
target_include_directories(Test
        PRIVATE "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Test PRIVATE Math)