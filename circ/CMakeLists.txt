set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_library(Circ OBJECT ${SOURCE_DIR}/pqc.c)
target_include_directories(Circ
        PRIVATE ${CMAKE_SOURCE_DIR}/include
        PUBLIC  ${INCLUDE_DIR}
)
target_link_libraries(Circ
        PRIVATE State
)

install(FILES       ${INCLUDE_DIR}/pqc.h
        DESTINATION include/qlib
)
