#include "special.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <boost/math/tools/precision.hpp>

using namespace std;

namespace {
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

    void test_legendre(unsigned const l, double const x, double const result) {
        TEST(legendre,  result, l, x);
        TEST(legendref, static_cast<float>(result), l, static_cast<float>(x));
        TEST(legendrel, static_cast<long double>(result), l, static_cast<long double>(x));
    }

    void test_hypot(double dx, double dy, double dz, double result) {
        TEST(hypot, result, dx, dy, dz);
        TEST(hypot, static_cast<float>(result),
            static_cast<float>(dx), static_cast<float>(dy), static_cast<float>(dz));
        TEST(hypot, static_cast<long double>(result),
            static_cast<long double>(dx), static_cast<long double>(dy), static_cast<long double>(dz));
    }
}

int main() {
    try {
        test_hypot(0.0, 0.0, 0.0, 0.0);
        test_hypot(1.0, 0.0, 0.0, 1.0);
        test_hypot(0.0, 1.0, 0.0, 1.0);
        test_hypot(0.0, 0.0, 1.0, 1.0);

        test_assoc_laguerre(0, 0, 0.0, 1.0);
        test_assoc_laguerre(0, 0, 0.5, 1.0);
        test_assoc_laguerre(0, 0, -0.5, 1.0);
        test_assoc_laguerre(0, 0, 1.0, 1.0);
        test_assoc_laguerre(0, 0, -1.0, 1.0);
        test_assoc_laguerre(0, 0, qNaN, qNaN);

        test_assoc_laguerre(0, 1, 0.0, 1.0);
        test_assoc_laguerre(0, 1, 0.5, 1.0);
        test_assoc_laguerre(0, 1, -0.5, 1.0);
        test_assoc_laguerre(0, 1, 1.0, 1.0);
        test_assoc_laguerre(0, 1, -1.0, 1.0);
        test_assoc_laguerre(0, 1, qNaN, qNaN);

        test_assoc_laguerre(1, 0, 0.0, 1.0);
        test_assoc_laguerre(1, 0, 0.5, 0.5);
        test_assoc_laguerre(1, 0, -0.5, 1.5);
        test_assoc_laguerre(1, 0, 1.0, 0.0);
        test_assoc_laguerre(1, 0, -1.0, 2.0);
        test_assoc_laguerre(1, 0, qNaN, qNaN);

        test_assoc_laguerre(1, 1, 0.0, 2.0);
        test_assoc_laguerre(1, 1, 0.5, 1.5);
        test_assoc_laguerre(1, 1, -0.5, 2.5);
        test_assoc_laguerre(1, 1, 1.0, 1.0);
        test_assoc_laguerre(1, 1, -1.0, 3.0);
        test_assoc_laguerre(1, 1, qNaN, qNaN);

        test_assoc_laguerre(2, 0, 0.0, 1.0);
        test_assoc_laguerre(2, 0, 0.5, 0.125);
        test_assoc_laguerre(2, 0, -0.5, 2.125);
        test_assoc_laguerre(2, 0, 1.0, -0.5);
        test_assoc_laguerre(2, 0, -1.0, 3.5);
        test_assoc_laguerre(2, 0, qNaN, qNaN);

        test_assoc_legendre(0, 0, 0.0, 1.0);
        test_assoc_legendre(0, 0, 0.5, 1.0);
        test_assoc_legendre(0, 0, -0.5, 1.0);
        test_assoc_legendre(0, 0, 1.0, 1.0);
        test_assoc_legendre(0, 0, -1.0, 1.0);
        test_assoc_legendre(0, 0, qNaN, qNaN);

        test_assoc_legendre(1, 0, 0.0, 0.0);
        test_assoc_legendre(1, 0, 0.5, 0.5);
        test_assoc_legendre(1, 0, -0.5, -0.5);
        test_assoc_legendre(1, 0, 1.0, 1.0);
        test_assoc_legendre(1, 0, -1.0, -1.0);
        test_assoc_legendre(1, 0, qNaN, qNaN);

        test_assoc_legendre(1, 1, 0.0, -1.0);
        test_assoc_legendre(1, 1, 0.5, -sqrt(3) / 2);
        test_assoc_legendre(1, 1, -0.5, -sqrt(3) / 2);
        test_assoc_legendre(1, 1, 1.0, 0.0);
        test_assoc_legendre(1, 1, -1.0, 0.0);
        test_assoc_legendre(1, 1, qNaN, qNaN);

        test_assoc_legendre(2, 0, 0.0, -0.5);
        test_assoc_legendre(2, 0, 0.5, -0.125);
        test_assoc_legendre(2, 0, -0.5, -0.125);
        test_assoc_legendre(2, 0, 1.0, 1.0);
        test_assoc_legendre(2, 0, -1.0, 1.0);
        test_assoc_legendre(2, 0, qNaN, qNaN);

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

        test_legendre(0, 0.0, 1.0);
        test_legendre(0, 0.5, 1.0);
        test_legendre(0, -0.5, 1.0);
        test_legendre(0, 1.0, 1.0);
        test_legendre(0, -1.0, 1.0);
        test_legendre(0, qNaN, qNaN);

        test_legendre(1, 0.0, 0.0);
        test_legendre(1, 0.5, 0.5);
        test_legendre(1, -0.5, -0.5);
        test_legendre(1, 1.0, 1.0);
        test_legendre(1, -1.0, -1.0);
        test_legendre(1, qNaN, qNaN);

        test_legendre(2, 0.0, -0.5);
        test_legendre(2, 0.5, -0.125);
        test_legendre(2, -0.5, -0.125);
        test_legendre(2, 1.0, 1.0);
        test_legendre(2, -1.0, 1.0);
        test_legendre(2, qNaN, qNaN);

        try {
            test_beta(0.0, 0.0, 0.0);
            std::cerr << "Expected exception\n";
            status = 1;
        } catch(std::domain_error const&) {}
    } catch(...) {
        std::cerr << "Accepted unexpected exception\n";
        status = 1;
    }

    return status;
}
