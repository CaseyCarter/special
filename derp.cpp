#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <random>
#include <boost/math/tools/precision.hpp>
#include "special.hpp"

namespace {
    template<class T>
    constexpr auto eps = boost::math::tools::epsilon<T>();

    template<class T>
    constexpr auto qNaN = std::numeric_limits<T>::quiet_NaN();

    template<class T>
    constexpr auto inf = std::numeric_limits<T>::infinity();

    auto const n = 500;
    auto engine = std::mt19937{};

    template<class T, T Deleter>
    struct deleter {
        template<class T>
        void operator()(T* ptr) noexcept {
            Deleter(ptr);
        }
    };

    using unique_file = std::unique_ptr<FILE, deleter<decltype(std::fclose), std::fclose>>;

    void generate(double const low, char const * const filename) {
        unique_file uf{std::fopen(filename, "wt")};
        FILE * const f = uf.get();
        assert(f);

        double const high = low / eps<double>;
        // std::printf("low: %.16e, high: %.16e\n", low, high);

        auto dist = std::uniform_int_distribution<unsigned long long>{};
        auto gen = [&] {
            auto in_range = [=](double const d) {
                return low <= d && d < high;
            };
            double d;
            do {
                auto const ull = dist(engine);
                std::memcpy(&d, &ull, sizeof(ull));
            } while(!in_range(std::fabs(d)));
            return d;
        };

        for (auto i = 0; i < n; ++i) {
            double const d[3] = {gen(), gen(), gen()};
            //std::printf("SC_(%+.16e), SC_(%+.16e), SC_(%+.16e), %+.16e\n",
            //    d[0], d[1], d[2], std::hypot(d[0], d[1], d[2]) - std::hypot(d[0], std::hypot(d[1], d[2])));
            std::fprintf(f, "sqrt(%.20e^2 + %.20e^2 + %.20e^2)  %.20e\n",
                std::fabs(d[0]), std::fabs(d[1]), std::fabs(d[2]), std::hypot(d[0], d[1], d[2]));
        }
    }
}

int main(int argc, char *argv[]) {
    generate(1.0, "hypot_low");
    generate(1.0e64, "hypot_high");
}
