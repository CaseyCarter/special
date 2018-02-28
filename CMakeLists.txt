cmake_minimum_required(VERSION 3.0)

project(special CXX)

add_library(boost INTERFACE)
target_include_directories(boost INTERFACE
    algorithm/include array/include assert/include bind/include
    concept_check/include config/include container/include container_hash/include
    core/include detail/include exception/include function/include functional/include
    integer/include io/include iterator/include lexical_cast/include math/include
    move/include mpl/include numeric_conversion/include optional/include predef/include
    preprocessor/include range/include smart_ptr/include static_assert/include
    test/include throw_exception/include timer/include type_index/include
    type_traits/include utility/include)
target_compile_definitions(boost INTERFACE
    BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE
    BOOST_MATH_DOMAIN_ERROR_POLICY=errno_on_error
    BOOST_MATH_OVERFLOW_ERROR_POLICY=ignore_error)
target_compile_options(boost INTERFACE /fp:precise)

enable_testing()

add_library(smf special_math.cpp)
target_link_libraries(smf PRIVATE boost)
target_compile_options(smf PRIVATE /std:c++latest /WX)

function(add_smf_test TESTNAME SOURCES...)
  add_executable(${TESTNAME} ${SOURCES...})
  target_compile_options(${TESTNAME} PRIVATE /std:c++latest /WX)
  target_link_libraries(${TESTNAME} smf boost)
  add_test(${TESTNAME} ${TESTNAME})
endfunction(add_smf_test)

add_smf_test(test_special test_special.cpp)
add_smf_test(test_wip test_wip.cpp)

