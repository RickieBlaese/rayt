#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <concepts>
#include <functional>
#include <numbers>
#include <optional>
#include <random>
#include <ratio>
#include <string>
#include <variant>
#include <vector>

#include <cinttypes>
#include <cmath>
#include <cstdio>

#include <notcurses/notcurses.h>

void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out);

float optional_min(std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc);

template <typename T>
T get_random_real(T a, T b) {
    static std::random_device device{};
    static std::default_random_engine engine(device());
    return std::uniform_real_distribution(a, b)(engine);
}


#define ERR_EXIT(...) { \
    std::fprintf(stderr, "error: file " __FILE__ ":%i in %s: ", __LINE__, __func__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    std::exit(EXIT_FAILURE); \
}

#endif
