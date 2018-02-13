cmake_minimum_required(VERSION 3.0)

project(special CXX)

add_library(boost INTERFACE)
target_include_directories(boost INTERFACE boostmath)

enable_testing()

add_library(smf
    assoc_laguerre.cpp assoc_laguerref.cpp
    assoc_legendre.cpp assoc_legendref.cpp
    beta.cpp betaf.cpp
    hypot3.cpp hypot3f.cpp
    legendre.cpp legendref.cpp)
target_link_libraries(smf PRIVATE boost)
target_compile_options(smf PUBLIC /std:c++latest)
target_compile_definitions(smf PUBLIC BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE)
target_include_directories(smf PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_executable(spot test.cpp)
target_link_libraries(spot smf boost)
add_test(spot spot)

add_executable(test_beta test_beta.cpp)
target_link_libraries(test_beta smf boost)
add_test(test_beta test_beta)
