#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <concepts>
#include <functional>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <ratio>
#include <string>
#include <variant>
#include <vector>

#include <cinttypes>
#include <cmath>
#include <cstring>
#include <cstdio>

#include <notcurses/notcurses.h>

void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out);

double optional_min(const std::pair<std::optional<double>, std::optional<double>> &ptimes, struct notcurses *nc);

template <typename T>
T get_random_real(T a, T b) {
    static thread_local std::random_device device{};
    static thread_local std::default_random_engine engine(device());
    return std::uniform_real_distribution(a, b)(engine);
}

template <typename T>
T get_random_int(T a, T b) {
    static thread_local std::random_device device{};
    static thread_local std::default_random_engine engine(device());
    return std::uniform_int_distribution(a, b)(engine);
}

using guid_t = std::uint64_t;

static constexpr guid_t sys_guid = std::numeric_limits<guid_t>::min();

guid_t generate_guid();


#define ERR_EXIT(...) { \
    std::fprintf(stderr, "error: file " __FILE__ ":%i in %s(): ", __LINE__, __func__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    std::exit(EXIT_FAILURE); \
}

#endif
