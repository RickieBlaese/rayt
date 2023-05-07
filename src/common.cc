#include "common.h"

void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out) {
    if (count <= 1) {
        out.emplace_back(a, b);
        return;
    }
    const float per = static_cast<float>(b - a) / static_cast<float>(count);

    std::int32_t i = 0;
    float x = 0.0f;
    float y = per;
    for (; i < count - 1; i++) {
        out.emplace_back(static_cast<std::int32_t>(std::round(x)),
            static_cast<std::int32_t>(std::round(y)));
        x += per;
        y += per;
    }
    out.emplace_back(static_cast<std::int32_t>(std::round(x)),
        static_cast<std::int32_t>(std::round(y)));
}

float optional_min(std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc) {
    if (ptimes.first.has_value() && !ptimes.second.has_value()) {
        return ptimes.first.value();
    }
    if (ptimes.second.has_value() && !ptimes.first.has_value()) {
        return ptimes.second.value();
    }
    if (!ptimes.first.has_value() && !ptimes.second.has_value()) {
        notcurses_stop(nc);
        ERR_EXIT("pair passed with two empty optionals");
    }
    return std::min(ptimes.first.value(), ptimes.second.value());
}
