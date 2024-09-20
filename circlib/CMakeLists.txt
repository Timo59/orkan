add_library(Circ OBJECT
        diag.c
        pauli.c
        circlib.c
        $<TARGET_OBJECTS:Math>
        $<TARGET_OBJECTS:State>
        $<TARGET_OBJECTS:Gate>
)

target_include_directories(Circ
        PRIVATE   "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Circ PRIVATE
        Q_compiler_flags
        Math
        State
        Gate
)

install(FILES "circlib.h" "diag.h" "pauli.h" DESTINATION include)