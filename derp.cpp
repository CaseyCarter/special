#define _CRT_SECURE_NO_WARNINGS

#include <cassert>
#include <cstdio>
#include <limits>
#include <random>
#include "special.hpp"

static_assert(sizeof(unsigned long long) == sizeof(double));

int main(int argc, char *argv[]) {
    auto n = 8;
    if (argc > 1) {
        std::sscanf(argv[1], "%d", &n);
        assert(n >= 0);
    }

    auto engine = std::mt19937{};
#if 0
    auto gen = [&engine, dist = std::uniform_real_distribution<>{std::numeric_limits<double>::min(), std::numeric_limits<double>::max()}]() mutable {
        double d;
        do {
            d = dist(engine);
        } while(std::isnan(d) || std::isinf(d));
        return d;
    };
#else
    auto gen = [&engine, dist = std::uniform_int_distribution<unsigned long long>{}]() mutable {
        double d;
        do {
            auto const ull = dist(engine);
            std::memcpy(&d, &ull, sizeof(ull));
        } while(std::isnan(d) || std::isinf(d));
        return d;
    };
#endif
    for (int i = 0; i < n; ++i) {
        double const d[3] = {gen(), gen(), gen()};
        std::printf("SC_(%.16e), SC_(%.16e), SC_(%.16e), SC_(%.16e)\n", d[0], d[1], d[2], std::hypot(d[0], d[1], d[2]));
    }
}
