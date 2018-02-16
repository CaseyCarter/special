#include <cmath>
#include <iostream>
#include <limits>
#include <boost/array.hpp>
#include "special.hpp"

namespace {
    int status = 0;

    template<class>
    constexpr const char *suffix = "";
    template<>
    constexpr const char *suffix<long double> = "l";
    template<>
    constexpr const char *suffix<float> = "f";

    template<class T>
    constexpr T relative_error(T const actual, T const expected) {
        if (expected)
            return static_cast<T>(std::abs(actual / expected - 1));
        else
            return static_cast<T>(std::abs(actual * 1e15));
    }

    template<class T, class First, class... Args, size_t... Is>
    void log_failure_(T const expected, T const actual, char const * const name,
        std::index_sequence<Is...>, First&& first, Args&&... args)
    {
        std::cerr.precision(12);
        std::cerr << "Check failed: " << name << suffix<T> << '(' << first;
        int unused[] = {0, ((std::cerr << ", " << args), 0)...};
        (void)unused;
        std::cerr << "): Actual: " << actual << ", Expected: " << expected
            << ", relative error: " << relative_error(actual, expected) << "\n\n";
        ++status;
        std::quick_exit(1);
    }

    template<class T, class First, class... Args>
    inline void log_failure(T const expected, T const actual, char const * const name,
        First&& first, Args&&... args)
    {
        log_failure_(expected, actual, name, std::index_sequence_for<Args...>{},
            std::forward<First>(first), std::forward<Args>(args)...);
    }

    template<class T>
    struct table_type { using type = T; };

    template<class T, class F>
    void test_beta_single(F f, T const x, T const y, T const expected, double const tolerance) {
        auto const actual = f(x, y);
        if (relative_error(actual, expected) > tolerance) {
            log_failure(expected, actual, "beta", x, y);
        }
    }

    template<class T, class F>
    void test_beta(F f) {
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))
#include "math/test/beta_small_data.ipp"
#include "math/test/beta_med_data.ipp"
#include "math/test/beta_exp_data.ipp"
#undef SC_

        for(auto& datum : beta_small_data) {
            test_beta_single(f, datum[0], datum[1], datum[2], 1e-15);
        }
        for(auto& datum : beta_med_data) {
            test_beta_single(f, datum[0], datum[1], datum[2], 1e-13);
        }
        for(auto& datum : beta_exp_data) {
            test_beta_single(f, datum[0], datum[1], datum[2], 1e-14);
        }
    }
} // unnamed namespace

int main() {
    test_beta<double>([](auto... args){ return std::beta(args...); });
    test_beta<float>([](auto... args){ return std::betaf(args...); });
    test_beta<long double>([](auto... args){ return std::betal(args...); });

    return status;
}
