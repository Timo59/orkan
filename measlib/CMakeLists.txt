add_library(Meas OBJECT
        exact.c
        $<TARGET_OBJECTS:Math>
        $<TARGET_OBJECTS:State>
        $<TARGET_OBJECTS:Gate>
        $<TARGET_OBJECTS:Circ>
)

target_include_directories(Meas
        PRIVATE   "${CMAKE_SOURCE_DIR}/include"
        INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Meas PRIVATE
        Q_compiler_flags
        Math
        State
        Gate
        Circ
)

install(FILES "exact.h" DESTINATION include)