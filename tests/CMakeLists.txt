cmake_minimum_required(VERSION 3.22)
project(Handler)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")

# Handler library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. ${CMAKE_BINARY_DIR}/Deexcitation)

enable_testing()

find_package(GTest REQUIRED)
include(GoogleTest)

set(Tests
    TestCOLA.cpp
    TestHandler.cpp
)

add_executable(RunTests ${Tests})

target_link_libraries(RunTests Deexcitation FermiBreakUp COLA)
target_link_libraries(RunTests GTest::GTest GTest::Main)

# default config
FILE(CREATE_LINK ${PROJECT_SOURCE_DIR}/data/config.xml config.xml)

gtest_discover_tests(RunTests)
