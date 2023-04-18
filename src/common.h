#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <functional>
#include <numbers>
#include <string>
#include <variant>
#include <vector>

#include <cinttypes>
#include <cmath>
#include <cstdio>

#include <notcurses/notcurses.h>

void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out);


#define ERR_EXIT(...) { \
    std::fprintf(stderr, "error: file " __FILE__ ", func %s, line %i: ", __func__, __LINE__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    std::exit(EXIT_FAILURE); \
}

#endif
