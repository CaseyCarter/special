#include "special.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>

using namespace std;

namespace {
    int status = 0;

    [[noreturn]] void fail(int const lineno, const char* const expr) {
        cerr << __FILE__ << '(' << lineno << "): Check failed: \"" << expr << "\"\n";
        status = 1;
    }

#define CHECK(...) ((__VA_ARGS__) ? void() : fail(__LINE__, #__VA_ARGS__))

    template<class T, class U, class C = common_type_t<T, U>>
    constexpr T narrow_cast(U const u) noexcept {
        auto const result = static_cast<T>(u);
        CHECK(C(result) == C(u)); // FIXME: cascade.
        return result;
    }

#define TEST(FN, L, M, X, RESULT)                                                             \
    do {                                                                                      \
        if (isnan(FN(L, M, X)) != isnan(result) || !isnan(result) && FN(L, M, X) != RESULT) { \
            cerr << "Check failed: " #FN "(" << L << ", " << M << ", " << X << "):\n"         \
                "Actual: " << FN(L, M, X) << ", Expected: " << RESULT << "\n\n";              \
            ::status = 1;                                                                     \
        }                                                                                     \
    } while(0)

    void test_assoc_laguerre(unsigned const l, unsigned const m, double const x, double const result) {
        TEST(assoc_laguerre,  l, m, x, result);
        TEST(assoc_laguerref, l, m, static_cast<float>(x), static_cast<float>(result));
        TEST(assoc_laguerrel, l, m, static_cast<long double>(x), static_cast<long double>(result));
    }

    void test_assoc_legendre(unsigned const l, unsigned const m, double const x, double const result) {
        TEST(assoc_legendre,  l, m, x, result);
        TEST(assoc_legendref, l, m, static_cast<float>(x), static_cast<float>(result));
        TEST(assoc_legendrel, l, m, static_cast<long double>(x), static_cast<long double>(result));
    }

#undef TEST

#define TEST(FN, L, X, RESULT)                                                \
    do {                                                                      \
        if (isnan(FN(L, X)) != isnan(result) ||                               \
           !isnan(result) && FN(L, X) != RESULT) {                            \
            cerr << "Check failed: " #FN "(" << L << ", " << X << "):\n"      \
                "Actual: " << FN(L, X) << ", Expected: " << RESULT << "\n\n"; \
            ::status = 1;                                                     \
        }                                                                     \
    } while(0)

    void test_legendre(unsigned const l, double const x, double const result) {
        TEST(legendre,  l, x, result);
        TEST(legendref, l, static_cast<float>(x), static_cast<float>(result));
        TEST(legendrel, l, static_cast<long double>(x), static_cast<long double>(result));
    }

#undef TEST
}

int main() {
    test_assoc_laguerre(0, 0, 0.0, 1.0);
    test_assoc_laguerre(0, 0, 0.5, 1.0);
    test_assoc_laguerre(0, 0, -0.5, 1.0);
    test_assoc_laguerre(0, 0, 1.0, 1.0);
    test_assoc_laguerre(0, 0, -1.0, 1.0);
    test_assoc_laguerre(0, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_laguerre(0, 1, 0.0, 1.0);
    test_assoc_laguerre(0, 1, 0.5, 1.0);
    test_assoc_laguerre(0, 1, -0.5, 1.0);
    test_assoc_laguerre(0, 1, 1.0, 1.0);
    test_assoc_laguerre(0, 1, -1.0, 1.0);
    test_assoc_laguerre(0, 1, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_laguerre(1, 0, 0.0, 1.0);
    test_assoc_laguerre(1, 0, 0.5, 0.5);
    test_assoc_laguerre(1, 0, -0.5, 1.5);
    test_assoc_laguerre(1, 0, 1.0, 0.0);
    test_assoc_laguerre(1, 0, -1.0, 2.0);
    test_assoc_laguerre(1, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_laguerre(1, 1, 0.0, 2.0);
    test_assoc_laguerre(1, 1, 0.5, 1.5);
    test_assoc_laguerre(1, 1, -0.5, 2.5);
    test_assoc_laguerre(1, 1, 1.0, 1.0);
    test_assoc_laguerre(1, 1, -1.0, 3.0);
    test_assoc_laguerre(1, 1, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_laguerre(2, 0, 0.0, 1.0);
    test_assoc_laguerre(2, 0, 0.5, 0.125);
    test_assoc_laguerre(2, 0, -0.5, 2.125);
    test_assoc_laguerre(2, 0, 1.0, -0.5);
    test_assoc_laguerre(2, 0, -1.0, 3.5);
    test_assoc_laguerre(2, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_legendre(0, 0, 0.0, 1.0);
    test_assoc_legendre(0, 0, 0.5, 1.0);
    test_assoc_legendre(0, 0, -0.5, 1.0);
    test_assoc_legendre(0, 0, 1.0, 1.0);
    test_assoc_legendre(0, 0, -1.0, 1.0);
    test_assoc_legendre(0, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_legendre(1, 0, 0.0, 0.0);
    test_assoc_legendre(1, 0, 0.5, 0.5);
    test_assoc_legendre(1, 0, -0.5, -0.5);
    test_assoc_legendre(1, 0, 1.0, 1.0);
    test_assoc_legendre(1, 0, -1.0, -1.0);
    test_assoc_legendre(1, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_legendre(1, 1, 0.0, -1.0);
    test_assoc_legendre(1, 1, 0.5, -sqrt(3) / 2);
    test_assoc_legendre(1, 1, -0.5, -sqrt(3) / 2);
    test_assoc_legendre(1, 1, 1.0, 0.0);
    test_assoc_legendre(1, 1, -1.0, 0.0);
    test_assoc_legendre(1, 1, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_assoc_legendre(2, 0, 0.0, -0.5);
    test_assoc_legendre(2, 0, 0.5, -0.125);
    test_assoc_legendre(2, 0, -0.5, -0.125);
    test_assoc_legendre(2, 0, 1.0, 1.0);
    test_assoc_legendre(2, 0, -1.0, 1.0);
    test_assoc_legendre(2, 0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_legendre(0, 0.0, 1.0);
    test_legendre(0, 0.5, 1.0);
    test_legendre(0, -0.5, 1.0);
    test_legendre(0, 1.0, 1.0);
    test_legendre(0, -1.0, 1.0);
    test_legendre(0, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_legendre(1, 0.0, 0.0);
    test_legendre(1, 0.5, 0.5);
    test_legendre(1, -0.5, -0.5);
    test_legendre(1, 1.0, 1.0);
    test_legendre(1, -1.0, -1.0);
    test_legendre(1, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());

    test_legendre(2, 0.0, -0.5);
    test_legendre(2, 0.5, -0.125);
    test_legendre(2, -0.5, -0.125);
    test_legendre(2, 1.0, 1.0);
    test_legendre(2, -1.0, 1.0);
    test_legendre(2, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
}
