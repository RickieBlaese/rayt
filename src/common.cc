#include "common.h"

void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out) {
    if (count <= 1) {
        out.emplace_back(a, b);
        return;
    }
    const float per = static_cast<float>(b - a) / static_cast<float>(count);

    std::int32_t i = 0;
    for (; i < count - 1; i++) {
        out.emplace_back(static_cast<std::int32_t>(std::round(static_cast<float>(i) * per + static_cast<float>(a))),
            static_cast<std::int32_t>(std::round(static_cast<float>(i + 1) * per + static_cast<float>(a))));
    }
    out.emplace_back(static_cast<float>(i) * per + static_cast<float>(a), b);
}
