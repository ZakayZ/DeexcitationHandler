cmake_minimum_required(VERSION 3.16)
project(Deexcitation VERSION 1.0.0 LANGUAGES CXX)

find_package(Geant4 REQUIRED)
find_package(COLA REQUIRED)

include(FetchContent)
FetchContent_Declare(
    FermiBreakUp
    GIT_REPOSITORY https://github.com/ZakayZ/FermiBreakUp.git
    GIT_TAG master
)
FetchContent_MakeAvailable(FermiBreakUp)

file(GLOB SOURCES
    Deexcitation/converter/*.cpp
    Deexcitation/factory/*.cpp
    Deexcitation/handler/*.cpp
)

add_library(Deexcitation ${SOURCES})

target_link_libraries(Deexcitation FermiBreakUp ${Geant4_LIBRARIES} COLA)

target_include_directories(Deexcitation PUBLIC ${Geant4_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})

target_compile_options(Deexcitation PRIVATE -Wall -Werror -Wfloat-conversion -Wextra -Wpedantic)

set_property(TARGET Deexcitation PROPERTY CXX_STANDARD 20)
set_property(TARGET Deexcitation PROPERTY CXX_STANDARD_REQUIRED On)

target_include_directories(Deexcitation PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:>
)

# Install the library and headers
install(TARGETS Deexcitation
    EXPORT DeexcitationExport
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(
    DIRECTORY Deexcitation
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
)

# Generate and install the package configuration file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DeexcitationConfig.cmake
    INSTALL_DESTINATION lib/cmake/Deexcitation
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/DeexcitationConfigVersion.cmake
    VERSION ${Deexcitation_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/DeexcitationConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/DeexcitationConfigVersion.cmake
    DESTINATION lib/cmake/Deexcitation
)

# Export the targets
install(EXPORT DeexcitationExport
    FILE DeexcitationExport.cmake
    DESTINATION lib/cmake/Deexcitation
)
