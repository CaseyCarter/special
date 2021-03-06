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
target_compile_options(boost INTERFACE /fp:precise /EHsc)

enable_testing()

add_library(smf special_math.cpp)
target_link_libraries(smf PRIVATE boost)
target_compile_options(smf PRIVATE /std:c++latest /WX)

add_executable(test_special test_special.cpp)
target_compile_options(test_special PRIVATE /std:c++latest /WX)
target_link_libraries(test_special smf boost)
add_test(test_special test_special)

# add_library(boost_multi INTERFACE)
# target_include_directories(boost_multi INTERFACE multiprecision/include rational/include)
# target_link_libraries(boost_multi INTERFACE boost)
#
# add_executable(gen_hypot_data gen_hypot_data.cpp)
# target_compile_options(gen_hypot_data PRIVATE /std:c++latest /WX)
# target_compile_definitions(gen_hypot_data PRIVATE _CRT_SECURE_NO_WARNINGS)
# target_link_libraries(gen_hypot_data boost_multi)
