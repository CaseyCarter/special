cmake_minimum_required(VERSION 3.0)

project(special CXX)

find_package(Boost REQUIRED)

enable_testing()

add_library(smf special.cpp)
target_include_directories(smf PRIVATE ${Boost_INCLUDE_DIR})
target_include_directories(smf PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_executable(spot test.cpp special.cpp)
target_include_directories(spot PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(spot smf)
add_test(spot spot)

add_executable(test_beta test_beta.cpp)
target_include_directories(test_beta PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(test_beta smf)
add_test(test_beta test_beta)
