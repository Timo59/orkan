# PlatformConfig.cmake - Platform-specific configuration

# Check the OS for specific build instructions
if(APPLE)
    message(STATUS "Building on macOS")
    # Prepend the libraries built with @rpath to set the runtime path for dyld
    set(CMAKE_MACOSX_RPATH 1)

elseif(UNIX)
    message(STATUS "Building on Linux")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

else()
    message(FATAL_ERROR "Unknown system: ${CMAKE_SYSTEM_NAME}")
endif()