cmake_minimum_required(VERSION 3.0)

project(special CXX)

find_package(Boost REQUIRED)

enable_testing()

add_executable(special test.cpp special.cpp)
target_include_directories(special PRIVATE ${Boost_INCLUDE_DIR})
add_test(special special)
