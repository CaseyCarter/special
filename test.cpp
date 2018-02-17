#define _USE_MATH_DEFINES

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <boost/math/tools/precision.hpp>
#include "special.hpp"

using namespace std;

namespace {
    constexpr auto inf = numeric_limits<double>::infinity();
    constexpr auto qNaN = numeric_limits<double>::quiet_NaN();

    int status = 0;

    template<class T, class First, class... Args, size_t... Is>
    void log_failure(T const expected, T const actual, char const * const name,
        std::index_sequence<Is...>, First&& first, Args&&... args)
    {
        cerr.precision(12);
        cerr << "Check failed: " << name << '(' << first;
        int unused[] = {0, ((cerr << ", " << args), 0)...};
        (void)unused;
        cerr << "): Actual: " << actual << ", Expected: " << expected << "\n\n";
        status = 1;
    }

    template<class T, class... Args>
    void test(T const expected, T const actual, char const * const name, Args... args) {
        if (isinf(expected)) {
            if (isinf(actual)) {
                return;
            }
        }
        if (isnan(expected)) {
            if (isnan(actual)) {
                return;
            }
        } else if (!isnan(actual)) {
            auto const absolute_error = expected - actual;
            auto const relative_error = abs(absolute_error / (expected ? expected : 1.0));
            if (relative_error < 1e-6) {
                return;
            }
        }
        log_failure(expected, actual, name, std::index_sequence_for<Args...>{}, args...);
    }

#define TEST(FN, EXPECTED, ...) test(EXPECTED, FN(__VA_ARGS__), #FN, __VA_ARGS__)

    void test_assoc_laguerre(unsigned const l, unsigned const m, double const x, double const result) {
        TEST(assoc_laguerre,  result, l, m, x);
        TEST(assoc_laguerref, static_cast<float>(result), l, m, static_cast<float>(x));
        TEST(assoc_laguerrel, static_cast<long double>(result), l, m, static_cast<long double>(x));
    }

    void test_assoc_legendre(unsigned const l, unsigned const m, double const x, double const result) {
        TEST(assoc_legendre,  result,  l, m, x);
        TEST(assoc_legendref, static_cast<float>(result), l, m, static_cast<float>(x));
        TEST(assoc_legendrel, static_cast<long double>(result), l, m, static_cast<long double>(x));
    }

    void test_beta(double const x, double const y, double const result) {
        TEST(beta,  result, x, y);
        TEST(betaf, static_cast<float>(result), static_cast<float>(x), static_cast<float>(y));
        TEST(betal, static_cast<long double>(result), static_cast<long double>(x), static_cast<long double>(y));

        TEST(beta,  result, y, x);
        TEST(betaf, static_cast<float>(result), static_cast<float>(y), static_cast<float>(x));
        TEST(betal, static_cast<long double>(result), static_cast<long double>(y), static_cast<long double>(x));
    }

    void test_comp_ellint_1(double const x, double const result) {
        TEST(comp_ellint_1,  result, x);
        TEST(comp_ellint_1f, static_cast<float>(result), static_cast<float>(x));
        TEST(comp_ellint_1l, static_cast<long double>(result), static_cast<long double>(x));
    }

    void test_comp_ellint_2(double const x, double const result) {
        TEST(comp_ellint_2,  result, x);
        TEST(comp_ellint_2f, static_cast<float>(result), static_cast<float>(x));
        TEST(comp_ellint_2l, static_cast<long double>(result), static_cast<long double>(x));
    }

    void test_comp_ellint_3(double const k, double const nu, double const result) {
        TEST(comp_ellint_3,  result, k, nu);
        TEST(comp_ellint_3f, static_cast<float>(result), static_cast<float>(k), static_cast<float>(nu));
        TEST(comp_ellint_3l, static_cast<long double>(result), static_cast<long double>(k), static_cast<long double>(nu));
    }

    void test_ellint_1(double const k, double const phi, double const result) {
        TEST(ellint_1,  result, k, phi);
        TEST(ellint_1f, static_cast<float>(result), static_cast<float>(k), static_cast<float>(phi));
        TEST(ellint_1l, static_cast<long double>(result), static_cast<long double>(k), static_cast<float>(phi));
    }

    void test_ellint_2(double const k, double const phi, double const result) {
        TEST(ellint_2,  result, k, phi);
        TEST(ellint_2f, static_cast<float>(result), static_cast<float>(k), static_cast<float>(phi));
        TEST(ellint_2l, static_cast<long double>(result), static_cast<long double>(k), static_cast<float>(phi));
    }

    void test_ellint_3(double const k, double const nu, double const phi, double const result) {
        TEST(ellint_3,  result, k, nu, phi);
        TEST(ellint_3f, static_cast<float>(result), static_cast<float>(k), static_cast<float>(nu), static_cast<float>(phi));
        TEST(ellint_3l, static_cast<long double>(result), static_cast<long double>(k), static_cast<float>(nu), static_cast<float>(phi));
    }

    void test_legendre(unsigned const l, double const x, double const result) {
        TEST(legendre,  result, l, x);
        TEST(legendref, static_cast<float>(result), l, static_cast<float>(x));
        TEST(legendrel, static_cast<long double>(result), l, static_cast<long double>(x));
    }

    void test_hypot(double dx, double dy, double dz, double result) {
        TEST(hypot, result, dx, dy, dz);
        TEST(hypot, result, dx, dz, dy);
        TEST(hypot, result, dz, dx, dy);
        TEST(hypot, result, dz, dy, dx);
        TEST(hypot, result, dy, dz, dx);
        TEST(hypot, result, dy, dx, dz);
        TEST(hypot, static_cast<float>(result),
            static_cast<float>(dx), static_cast<float>(dy), static_cast<float>(dz));
        TEST(hypot, static_cast<long double>(result),
            static_cast<long double>(dx), static_cast<long double>(dy), static_cast<long double>(dz));
    }
}

int main() {
    try {
        test_assoc_laguerre(0, 0,  0.0,  1.0);
        test_assoc_laguerre(0, 0,  0.5,  1.0);
        test_assoc_laguerre(0, 0, -0.5,  1.0);
        test_assoc_laguerre(0, 0,  1.0,  1.0);
        test_assoc_laguerre(0, 0, -1.0,  1.0);
        test_assoc_laguerre(0, 0, qNaN, qNaN);

        test_assoc_laguerre(0, 1,  0.0,  1.0);
        test_assoc_laguerre(0, 1,  0.5,  1.0);
        test_assoc_laguerre(0, 1, -0.5,  1.0);
        test_assoc_laguerre(0, 1,  1.0,  1.0);
        test_assoc_laguerre(0, 1, -1.0,  1.0);
        test_assoc_laguerre(0, 1, qNaN, qNaN);

        test_assoc_laguerre(1, 0,  0.0,  1.0);
        test_assoc_laguerre(1, 0,  0.5,  0.5);
        test_assoc_laguerre(1, 0, -0.5,  1.5);
        test_assoc_laguerre(1, 0,  1.0,  0.0);
        test_assoc_laguerre(1, 0, -1.0,  2.0);
        test_assoc_laguerre(1, 0, qNaN, qNaN);

        test_assoc_laguerre(1, 1,  0.0,  2.0);
        test_assoc_laguerre(1, 1,  0.5,  1.5);
        test_assoc_laguerre(1, 1, -0.5,  2.5);
        test_assoc_laguerre(1, 1,  1.0,  1.0);
        test_assoc_laguerre(1, 1, -1.0,  3.0);
        test_assoc_laguerre(1, 1, qNaN, qNaN);

        test_assoc_laguerre(2, 0,  0.0,   1.0);
        test_assoc_laguerre(2, 0,  0.5, 0.125);
        test_assoc_laguerre(2, 0, -0.5, 2.125);
        test_assoc_laguerre(2, 0,  1.0,  -0.5);
        test_assoc_laguerre(2, 0, -1.0,   3.5);
        test_assoc_laguerre(2, 0, qNaN,  qNaN);

        test_assoc_legendre(0, 0,  0.0,  1.0);
        test_assoc_legendre(0, 0,  0.5,  1.0);
        test_assoc_legendre(0, 0, -0.5,  1.0);
        test_assoc_legendre(0, 0,  1.0,  1.0);
        test_assoc_legendre(0, 0, -1.0,  1.0);
        test_assoc_legendre(0, 0, qNaN, qNaN);

        test_assoc_legendre(1, 0,  0.0,  0.0);
        test_assoc_legendre(1, 0,  0.5,  0.5);
        test_assoc_legendre(1, 0, -0.5, -0.5);
        test_assoc_legendre(1, 0,  1.0 , 1.0);
        test_assoc_legendre(1, 0, -1.0, -1.0);
        test_assoc_legendre(1, 0, qNaN, qNaN);

        test_assoc_legendre(1, 1,  0.0,         -1.0);
        test_assoc_legendre(1, 1,  0.5, -sqrt(3) / 2);
        test_assoc_legendre(1, 1, -0.5, -sqrt(3) / 2);
        test_assoc_legendre(1, 1,  1.0,          0.0);
        test_assoc_legendre(1, 1, -1.0,          0.0);
        test_assoc_legendre(1, 1, qNaN,         qNaN);

        test_assoc_legendre(2, 0,  0.0,   -0.5);
        test_assoc_legendre(2, 0,  0.5, -0.125);
        test_assoc_legendre(2, 0, -0.5, -0.125);
        test_assoc_legendre(2, 0,  1.0,    1.0);
        test_assoc_legendre(2, 0, -1.0,    1.0);
        test_assoc_legendre(2, 0, qNaN,   qNaN);

        test_assoc_legendre(2, 2,  0.0,  3.0);
        test_assoc_legendre(2, 2,  0.5, 2.25);
        test_assoc_legendre(2, 2, -0.5, 2.25);
        test_assoc_legendre(2, 2,  1.0,  0.0);
        test_assoc_legendre(2, 2, -1.0,  0.0);
        test_assoc_legendre(2, 2, qNaN, qNaN);

        test_beta(1.0, 1.0, 1.0/1);
        test_beta(1.0, 2.0, 1.0/2);
        test_beta(1.0, 3.0, 1.0/3);
        test_beta(1.0, 4.0, 1.0/4);
        test_beta(1.0, 5.0, 1.0/5);

        test_beta(2.0, 2.0, 1.0/6);
        test_beta(2.0, 3.0, 1.0/12);
        test_beta(2.0, 4.0, 1.0/20);
        test_beta(2.0, 5.0, 1.0/30);

        test_beta(3.0, 3.0, 1.0/30);
        test_beta(3.0, 4.0, 1.0/60);
        test_beta(3.0, 5.0, 1.0/105);

        test_beta(4.0, 4.0, 1.0/140);
        test_beta(4.0, 5.0, 1.0/280);

        test_beta(5.0, 5.0, 1.0/630);

        test_beta(1.0, qNaN, qNaN);
        test_beta(qNaN, 1.0, qNaN);

        test_beta(1.0, 4.0, 0.25);
        {
            const auto small = boost::math::tools::epsilon<double>() / 1024;
            test_beta(small, 4.0, 1/small);
        }
        test_beta(4.0, 20.0, 0.00002823263692828910220214568040654997176736);
        test_beta(0.0125, 0.000023, 43558.24045647538375006349016083320744662);
        try {
            test_beta(0.0, 0.0, 0.0);
            std::cerr << "Expected exception\n";
            status = 1;
        } catch(std::domain_error const&) {}

        test_comp_ellint_1( 0.0,  M_PI_2);
        test_comp_ellint_1( 0.5, 1.68575);
        test_comp_ellint_1(-0.5, 1.68575);
        test_comp_ellint_1(qNaN,    qNaN);

        test_comp_ellint_2( 0.0, M_PI_2);
        test_comp_ellint_2( 1.0,    1.0);
        test_comp_ellint_2(qNaN,   qNaN);

        test_comp_ellint_3( 0.0,  0.0,  M_PI_2);
        test_comp_ellint_3( 0.5,  0.0, 1.68575);
        test_comp_ellint_3(qNaN,  0.0,    qNaN);
        test_comp_ellint_3( 0.0, qNaN,    qNaN);

        test_ellint_1( 0.0,  M_PI_2,   M_PI_2);
        test_ellint_1( 0.0, -M_PI_2,  -M_PI_2);
        test_ellint_1( 0.5,  M_PI_2,  1.68575);
        test_ellint_1(-0.5,  M_PI_2,  1.68575);
        test_ellint_1(-0.5, -M_PI_2, -1.68575);
        test_ellint_1( 0.7,       0,        0);
        test_ellint_1(-0.7,       0,        0);
        test_ellint_1(qNaN,  M_PI_2,     qNaN);
        test_ellint_1( 0.0,    qNaN,     qNaN);

        test_ellint_2( 0.0,  M_PI_2,  M_PI_2);
        test_ellint_2( 0.0, -M_PI_2, -M_PI_2);
        test_ellint_2( 0.7,     0.0,     0.0);
        test_ellint_2(-0.7,     0.0,     0.0);
        test_ellint_2( 1.0,  M_PI_2,     1.0);
        test_ellint_2( 1.0, -M_PI_2,    -1.0);
        test_ellint_2(qNaN,  M_PI_2,    qNaN);
        test_ellint_2( 0.0,    qNaN,    qNaN);

        test_ellint_3( 0.0,  0.0, M_PI_2,  M_PI_2);
        test_ellint_3( 0.5,  0.0, M_PI_2, 1.68575);
        test_ellint_3(qNaN,  0.0, M_PI_2,    qNaN);
        test_ellint_3( 0.0, qNaN, M_PI_2,    qNaN);

        test_hypot(0.0,   0.0, 0.0, 0.0);
        test_hypot(1.0,   0.0, 0.0, 1.0);
        test_hypot(+inf,  0.0, 1.0, inf);
        test_hypot(-inf,  0.0, 1.0, inf);
        test_hypot(+inf, qNaN, 1.0, inf); // C11 F.10.4.3: "hypot(+/-inf, y)
        test_hypot(-inf, qNaN, 1.0, inf); // returns +inf even if y is NaN"

        test_legendre(0,  0.0,  1.0);
        test_legendre(0,  0.5,  1.0);
        test_legendre(0, -0.5,  1.0);
        test_legendre(0,  1.0,  1.0);
        test_legendre(0, -1.0,  1.0);
        test_legendre(0, qNaN, qNaN);

        test_legendre(1,  0.0,  0.0);
        test_legendre(1,  0.5,  0.5);
        test_legendre(1, -0.5, -0.5);
        test_legendre(1,  1.0,  1.0);
        test_legendre(1, -1.0, -1.0);
        test_legendre(1, qNaN, qNaN);

        test_legendre(2,  0.0,   -0.5);
        test_legendre(2,  0.5, -0.125);
        test_legendre(2, -0.5, -0.125);
        test_legendre(2,  1.0,    1.0);
        test_legendre(2, -1.0,    1.0);
        test_legendre(2, qNaN,   qNaN);
    } catch(...) {
        std::cerr << "Caught unexpected exception\n";
        status = 1;
    }

    return status;
}
