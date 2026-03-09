cmake_minimum_required(VERSION 3.16)
project(Deexcitation VERSION 1.0.0 LANGUAGES CXX)

find_package(Geant4 REQUIRED)
find_package(COLA REQUIRED)

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME Deexcitation)

include(FetchContent)
include(GNUInstallDirs)

file(GLOB SOURCES
    Deexcitation/converter/*.cpp
    Deexcitation/factory/*.cpp
    Deexcitation/handler/*.cpp
)

add_library(${COLA_MODULE_NAME} SHARED ${SOURCES})

target_link_libraries(${COLA_MODULE_NAME} COLA ${Geant4_LIBRARIES})
target_include_directories(${COLA_MODULE_NAME} PUBLIC ${Geant4_INCLUDE_DIR})

target_compile_options(${COLA_MODULE_NAME} PRIVATE -Wall -Werror -Wextra -Wpedantic)

target_compile_features(${COLA_MODULE_NAME} PRIVATE cxx_std_17)

target_include_directories(${COLA_MODULE_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Deexcitation>
)

# Install the library and headers
install(TARGETS ${COLA_MODULE_NAME}
    EXPORT ${COLA_MODULE_NAME}Export
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY ${COLA_MODULE_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.h"
)

# Generate and install the package configuration file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
    VERSION ${Deexcitation_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
)

# Export the targets
install(EXPORT ${COLA_MODULE_NAME}Export
    FILE ${COLA_MODULE_NAME}Export.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
)
