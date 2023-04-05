#include <iostream>
#include <utility>
#include <chrono>
#include <numbers>
#include <map>
#include <vector>
#include <algorithm>

#include <cinttypes>
#include <cmath>

#include <notcurses/notcurses.h>


#define ERR_EXIT(...) { \
    fprintf(stderr, "error: file " __FILE__ ", func %s, line %i: ", __func__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fputc('\n', stderr); \
    exit(EXIT_FAILURE); \
}


template <typename T>
struct gvec3_t { /* NOLINT */
    T x = 0, y = 0, z = 0;

    gvec3_t() = default;
    gvec3_t(T x, T y, T z) : x(x), y(y), z(z) {}
    gvec3_t(const gvec3_t& other) = default;
    ~gvec3_t() = default;
    bool operator==(const gvec3_t& other) const { return x == other.x && y == other.y && z == other.z; }
    gvec3_t& operator=(const gvec3_t& other) = default;

    gvec3_t operator+(const gvec3_t& other) const { return {x + other.x, y + other.y, z + other.z}; }
    gvec3_t& operator+=(const gvec3_t& other) { return *this = *this + other; }
    gvec3_t operator-() const { return {-x, -y, -z}; }
    gvec3_t operator-(const gvec3_t& other) const { return *this + (-other); }
    gvec3_t& operator-=(const gvec3_t& other) { return *this = *this - other; }
    template <typename Y>
    gvec3_t<T> operator*(const Y v) const { return {x * v, y * v, z * v}; }
    template <typename Y>
    gvec3_t<T> operator/(const Y v) const { return {x / v, y / v, z / v}; }
    float mod() const { return std::sqrt(x * x + y * y + z * z); }
    T dot(const gvec3_t& other) const { return x * other.x + y * other.y + z * other.z; }
    gvec3_t normalized() const { return *this / mod(); }
};

template <typename T>
struct gvec2_t { /* NOLINT */
    T x = 0, z = 0;

    gvec2_t() = default;
    gvec2_t(T x, T z) : x(x), z(z) {}
    gvec2_t(const gvec2_t& other) = default;
    ~gvec2_t() = default;
    bool operator==(const gvec2_t& other) const { return x == other.x && z == other.z; }
    gvec2_t& operator=(const gvec2_t& other) = default;

    gvec2_t operator+(const gvec2_t& other) const { return {x + other.x, z + other.z}; }
    gvec2_t& operator+=(const gvec2_t& other) { return *this = *this + other; }
    gvec2_t operator-() const { return {-x, -z}; }
    gvec2_t operator-(const gvec2_t& other) const { return *this + (-other); }
    gvec2_t& operator-=(const gvec2_t& other) { return *this = *this - other; }
    gvec2_t operator*(const T v) const { return {x * v, z * v}; }
    gvec2_t operator/(const T v) const { return {x / v, z / v}; }

    float mod() const { return std::sqrt(x * x + z * z); }
    T dot(const gvec2_t& other) const { return x * other.x + z * other.z; }
    gvec2_t normalized() const { return *this / mod(); }
};


using vec3_t = gvec3_t<float>;
using vec2_t = gvec2_t<float>;
using rgb_t = gvec3_t<std::int32_t>;

/*
vec3_t operator"" _vec(const char *s) {
    std::string in = s;
    if (in.find('{') == std::string::npos || in.find('}') == std::string::npos || (in.find('{') > in.find('}'))) {
        throw std::runtime_error("cannot convert \"" + in + "\" to vec3_t");
    }
    in = in.substr(in.find('{'), in.find('}') - in.find('}'));
    std::vector<std::int32_t> nums;
    for (decltype(in)::iterator it = in.begin(); it != in.end(); it++) {
        nums.push_back(
    return {std::atoi(in
} */



struct line_t {
    vec3_t pos;
    float a = 0, b = 0, c = 0;

    vec3_t f(float t) const {
        return {pos.x + a * t, pos.y + b * t, pos.z + c * t};
    }
};

line_t line_between(const vec3_t& p1, const vec3_t& p2) {
    return {p1, p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
}


/* not even sure if this is the correct representation */
struct plane_t {
    vec2_t pos;
    float a = 0, b = 0;

    vec3_t f(float t, float u) const {
        const float x = pos.x + a * t;
        const float z = pos.z + b * u;
        const float y = std::sqrt(t * t - x * x) + std::sqrt(u * u - z * z);
        return {x, y, z};
    }
};

plane_t plane_between(const vec3_t& p1, const vec3_t& p2, const vec3_t& p3) {
    return {vec2_t{p1.x, p1.z}, p2.x - p1.x, p3.y - p2.y};
}


struct sphere_t {
    vec3_t pos;
    float r = 0;
};

void ls_intersect(const line_t& line, const sphere_t& sphere, std::vector<vec3_t>& out) {
    const float alpha = (line.a * line.a) + (line.b * line.b) + (line.c * line.c);
    const float beta  = 2 * (
        line.a * (line.pos.x - sphere.pos.x) +
        line.b * (line.pos.y - sphere.pos.y) +
        line.c * (line.pos.z - sphere.pos.z));
    const float gamma = -2 * ((line.pos.x * sphere.pos.x) + (line.pos.y * sphere.pos.y) + (line.pos.z * sphere.pos.z)) +
        line.pos.x * line.pos.x + sphere.pos.x * sphere.pos.x +
        line.pos.y * line.pos.y + sphere.pos.y * sphere.pos.y +
        line.pos.z * line.pos.z + sphere.pos.z * sphere.pos.z -
        (sphere.r * sphere.r);
    const float discr = (beta * beta) - 4 * alpha * gamma;
    if (discr < 0) {
        return;
    }
    const float t1    = (-beta + std::sqrt(discr)) / (2 * alpha);
    const float t2    = (-beta - std::sqrt(discr)) / (2 * alpha);
    out.push_back(line.f(t1));
    out.push_back(line.f(t2));
}

inline void sort_by_dist(std::vector<vec3_t>& vecs, const vec3_t& pos) {
    std::sort(vecs.begin(), vecs.end(), [&pos](const vec3_t& a, const vec3_t& b) { return (pos - a).mod() < (pos - b).mod(); });
}

struct vec4_t {
    vec3_t a, b;
};

/* acts like point light */
struct light_t {
    sphere_t sphere;
    float strength = 0.0f;
};

struct obj_t {
    sphere_t sphere;
    float smoothness = 1.0f;
    rgb_t color;
};

struct scene_t {
    std::vector<obj_t> objects;
    std::vector<light_t> lights;
};

struct char_ex_info_t { /* NOLINT */
    char ch = ' ';
    rgb_t fg{255, 255, 255}, bg{22, 24, 33};

    char_ex_info_t& operator=(const char_ex_info_t& other) = default;

    bool operator==(const char_ex_info_t& other) const {
        return ch == other.ch && fg == other.fg && bg == other.bg;
    }

    bool operator!=(const char_ex_info_t& other) const {
        return !(*this == other);
    }
};

std::string move_out(std::int32_t y, std::int32_t x) {
    return "\033[" + std::to_string(y) + ';' + std::to_string(x) + 'H';
}


std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}


int main() {
    if (!setlocale(LC_ALL, "")) {
        ERR_EXIT("couldn't set locale");
    }
    constexpr float fps = 60.0f;

    const std::string gradient = ".._,'`^\"-~:;=l!i><+?|][}{)(\\/1trxnuvczjfXYUJICLQO0Zmwqpdbkhao*#MW&8%B@$";
    const auto gradient_length = static_cast<std::int32_t>(gradient.length());
    const auto get_gradient = [&gradient, &gradient_length](float x) -> char {
        return gradient[std::clamp<std::int32_t>(static_cast<std::int32_t>(std::round(x)), 0, gradient_length - 1)];
    };

    const auto wait_us_per_frame = static_cast<std::uint64_t>(1'000'000.0 / fps);

    /* notcurses initialize */
    notcurses_options opts{};
    opts.flags = NCOPTION_INHIBIT_SETLOCALE
        | NCOPTION_DRAIN_INPUT;

    struct notcurses *nc = notcurses_init(nullptr, nullptr);

    if (nc == nullptr) {
        ERR_EXIT("couldn't initialize notcurses, notcurses_init returned nullptr");
    }

    struct ncplane *plane = notcurses_stdplane(nc);
    ncplane_erase(plane);

    vec3_t camera_pos(0, 0, 0);

    /* gaussian distribution */
    /* const auto gaussian = [](float x, float y) { return std::pow(std::numbers::e, -x * x - y * y); }; */

    std::uint64_t begin_time = get_current_time();
    std::uint64_t last_time  = begin_time;
    const float begin_draw_dist = 30.0f;

    /* aspect ratio ? */
    const float x_mul = 1.0f;
    const float y_mul = 1.0f;

    const float z_move_speed = 1.0f;
    const float x_move_speed = 1.0f;
    const float y_move_speed = 1.0f;

    scene_t scene;
    scene.objects.push_back(obj_t{
        sphere_t{vec3_t(0, 0, 30), 10.0f},
        0.01f,
        {255, 255, 255}
    });
    scene.lights.push_back(light_t{
        sphere_t{vec3_t(0, 15, 30), 1.0f},
        0.5f
    });
    const vec3_t light_orig_pos = scene.lights[0].sphere.pos;

    std::vector<vec3_t> ipoints;
    std::map<float, char_ex_info_t> dist_to_chars;
    std::vector<std::vector<char_ex_info_t>> buffer, last_buffer;

    std::uint32_t dimy = 0, dimx = 0;
    ncplane_dim_yx(plane, &dimy, &dimx);
    buffer.resize(dimy);
    for (std::int32_t i = 0; i < dimy; i++) {
        buffer[i].resize(dimx);
    }
    last_buffer = buffer;

    std::vector<std::tuple<std::int32_t, std::int32_t, char>> updated_positions;

    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        last_time = get_current_time();

        ncplane_dim_yx(plane, &dimy, &dimx);


        /* resize buffers as window size changes */
        if (buffer.size() != dimy) {
            buffer.resize(dimy);
        }

        if (buffer.size() > 0) {
            if (buffer[0].size() != dimx) {
                for (std::vector<char_ex_info_t>& line : buffer) {
                    line.resize(dimx);
                }
            }
        }

        /* revolve the light source */
        scene.lights[0].sphere.pos.x = light_orig_pos.x + 10 * std::cos(static_cast<float>(last_time - begin_time) / 1'000'000.0f);
        scene.lights[0].sphere.pos.z = light_orig_pos.z + 10 * std::sin(static_cast<float>(last_time - begin_time) / 1'000'000.0f);

        
        uint32_t chin = notcurses_get_nblock(nc, nullptr);
        if (chin == L'\t') { break; }
        else if (chin == L'q') { camera_pos.y += y_move_speed; }
        else if (chin == L'e') { camera_pos.y -= y_move_speed; }
        else if (chin == L'w') { camera_pos.z += z_move_speed; }
        else if (chin == L's') { camera_pos.z -= z_move_speed; }
        else if (chin == L'd') { camera_pos.x += x_move_speed; }
        else if (chin == L'a') { camera_pos.x -= x_move_speed; }


        for (std::int32_t i = 0; i < dimy; i++) {
            for (std::int32_t j = 0; j < dimx / 2; j++) {
                /* position for i is flipped since ncurses says y down is positive while we want y up is positive */
                vec3_t curpos = vec3_t((j - dimx / 4.0f) * x_mul, (-i + dimy / 2.0f) * y_mul, begin_draw_dist); /* NOLINT */
                line_t ray = line_between(camera_pos, camera_pos + curpos);
                dist_to_chars.clear();

                /* render spheres with appropriate lighting */
                for (const obj_t& object : scene.objects) {
                    ipoints.clear();
                    ls_intersect(ray, object.sphere, ipoints);
                    if (!ipoints.empty()) {
                        sort_by_dist(ipoints, camera_pos);
                        vec3_t sphere_minvec = ipoints[0];
                        vec3_t sphere_maxvec = ipoints[1];
                        /* nonpermanent and bad solution for camera being in front of sphere */
                        if (sphere_minvec.z > camera_pos.z) {
                            /* checking to make sure it is the visible ls_intersection
                             * i.e. the closest ls_intersect with sphere is equal or close
                             * enough to the closest ls_intersect with light source */
                            float applied_light = 0.0f;
                            for (const light_t& light : scene.lights) {
                                ipoints.clear();
                                const float dotp = (object.sphere.pos - sphere_minvec).normalized().dot((object.sphere.pos - light.sphere.pos).normalized());
                                ls_intersect(line_between(light.sphere.pos, sphere_minvec), object.sphere, ipoints);
                                /* guaranteed to ls_intersect because by above line it goes through both light and sphere */
                                vec3_t lminvec;
                                if (ipoints.empty()) {
                                    lminvec = sphere_minvec;
                                } else {
                                    sort_by_dist(ipoints, light.sphere.pos);
                                    lminvec = ipoints[0];
                                }
                                if (dotp > object.smoothness && (lminvec - sphere_minvec).mod() < 0.05f) {
                                    applied_light += light.strength * dotp;
                                }
                            }
                            applied_light = std::clamp<float>(applied_light, 0.0f, 1.0f);
                            char outch = '.';
                            if (applied_light > 0.0f) {
                                outch = get_gradient(static_cast<float>(gradient_length - 1) * applied_light);
                            } else if ((sphere_minvec - sphere_maxvec).mod() < 0.1f) { /* is close, i.e. edge */
                                outch = '-';
                            }
                            gvec3_t<std::int32_t> outcolor;
                            outcolor.x = static_cast<std::int32_t>(static_cast<float>(object.color.x) * applied_light);
                            outcolor.y = static_cast<std::int32_t>(static_cast<float>(object.color.y) * applied_light);
                            outcolor.z = static_cast<std::int32_t>(static_cast<float>(object.color.z) * applied_light);
                            dist_to_chars[(camera_pos - sphere_minvec).mod()] = char_ex_info_t{outch, outcolor};
                        }
                    }
                }

                /* render light sources */
                for (const light_t& light : scene.lights) {
                    /* checking ls_intersect with light src */
                    ipoints.clear();
                    ls_intersect(ray, light.sphere, ipoints);
                    if (!ipoints.empty()) {
                        /* find closest point of ls_intersect, we do not want back of sphere */
                        sort_by_dist(ipoints, camera_pos);
                        const vec3_t minvec = ipoints[0];

                        /* nonpermanent and bad solution */
                        if (minvec.z > camera_pos.z) {
                            dist_to_chars[(camera_pos - light.sphere.pos).mod()] = char_ex_info_t{get_gradient(light.strength)};
                        }
                    }
                }

                if (dist_to_chars.empty()) {
                    buffer[i][j] = char_ex_info_t{' '};
                } else {
                    buffer[i][j] = dist_to_chars.begin()->second;
                }
            }
        }

        for (std::int32_t i = 0; i < buffer.size(); i++) {
            for (std::int32_t j = 0; j < buffer[i].size(); j++) {
                /* these ifs are separated so as not to rely on short circuiting */
                if (i < last_buffer.size()) {
                    if (j < last_buffer[i].size()) {
                        if (buffer[i][j] != last_buffer[i][j]) {
                            updated_positions.emplace_back(i, j, buffer[i][j].ch);
                            if (updated_positions.size() > dimy - 3) {
                                updated_positions.erase(updated_positions.begin() + static_cast<std::int64_t>((updated_positions.size() - (dimy - 3))));
                            }
                            ncplane_set_fg_rgb8(plane, buffer[i][j].fg.x, buffer[i][j].fg.y, buffer[i][j].fg.z);
                            /* ncplane_set_bg_rgb8(plane, buffer[i][j].bg.x, buffer[i][j].bg.y, buffer[i][j].bg.z); */
                            ncplane_putwc_yx(plane, i, j * 2, buffer[i][j].ch);
                            ncplane_putwc(plane, buffer[i][j].ch);
                        }
                    }
                }
            }
        }
        last_buffer = buffer;

        /* mvwprintw(win, 0, 0, "x: %.2f", camera_pos.x);
        mvwprintw(win, 1, 0, "y: %.2f", camera_pos.y);
        mvwprintw(win, 2, 0, "z: %.2f", camera_pos.z);
        std::int32_t cury = 2;
        for (const auto& [x, y, ch] : updated_positions) {
            if (cury >= LINES) { break; }
            cury++;
            mvwprintw(win, cury, 0, "%c, x: %4i, y: %-4i", ch, x, y);
        }

        wrefresh(win); */
        notcurses_render(nc);
    }


    ncplane_destroy(plane);
    notcurses_stop(nc);
    return 0;
}
