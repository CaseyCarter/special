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
target_compile_definitions(boost INTERFACE BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE)

enable_testing()

add_library(smf
    assoc_laguerre.cpp assoc_laguerref.cpp
    assoc_legendre.cpp assoc_legendref.cpp
    beta.cpp betaf.cpp
    comp_ellint1.cpp comp_ellint1f.cpp
    comp_ellint2.cpp comp_ellint2f.cpp
    comp_ellint3.cpp comp_ellint3f.cpp
    cyl_bessel_i.cpp cyl_bessel_if.cpp
    cyl_bessel_j.cpp cyl_bessel_jf.cpp
    cyl_bessel_k.cpp cyl_bessel_kf.cpp
    cyl_neumann.cpp cyl_neumannf.cpp
    ellint1.cpp ellint1f.cpp
    ellint2.cpp ellint2f.cpp
    ellint3.cpp ellint3f.cpp
    expint.cpp expintf.cpp
    hermite.cpp hermitef.cpp
    hypot3.cpp hypot3f.cpp
    laguerre.cpp laguerref.cpp
    legendre.cpp legendref.cpp
    riemann_zeta.cpp riemann_zetaf.cpp
    sph_bessel.cpp sph_besself.cpp
    sph_legendre.cpp sph_legendref.cpp
    sph_neumann.cpp sph_neumannf.cpp)
target_link_libraries(smf PRIVATE boost)
target_compile_options(smf PRIVATE /std:c++latest)

add_executable(spot test.cpp)
target_compile_options(spot PRIVATE /std:c++latest)
target_link_libraries(spot smf boost)
add_test(spot spot)

add_executable(test_beta test_beta.cpp)
target_compile_options(test_beta PRIVATE /std:c++latest)
target_link_libraries(test_beta smf boost)
add_test(test_beta test_beta)
