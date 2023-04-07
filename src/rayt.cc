#include <iostream>
#include <limits>
#include <unistd.h>
#include <utility>
#include <chrono>
#include <numbers>
#include <functional>
#include <map>
#include <unordered_map>
#include <variant>
#include <vector>
#include <random>
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


/* general vec3 */
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
    gvec3_t<T> operator*(const Y& v) const { return {x * v, y * v, z * v}; }
    template <typename Y>
    gvec3_t<T> operator/(const Y& v) const { return {x / v, y / v, z / v}; }

    float mod() const { return std::sqrt(x * x + y * y + z * z); }
    T dot(const gvec3_t& other) const { return x * other.x + y * other.y + z * other.z; }
    gvec3_t normalized() const { return *this / mod(); }
    gvec3_t x_rotated(float theta) const { return {x + mod() * std::cos(theta), y, z}; }
    gvec3_t y_rotated(float theta) const { return {x, y + mod() * std::sin(theta), z}; }
    gvec3_t z_rotated(float theta) const { return {x, y, z + mod() * std::sin(theta)}; }
};

/* general vec2 */
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


inline void sort_by_dist(std::vector<vec3_t>& vecs, const vec3_t& pos) {
    std::sort(vecs.begin(), vecs.end(), [&pos](const vec3_t& a, const vec3_t& b) { return (pos - a).mod() < (pos - b).mod(); });
}


struct line_t {
    vec3_t pos, n;

    vec3_t f(float t) const {
        return {pos.x + n.x * t, pos.y + n.y * t, pos.z + n.z * t};
    }
};

line_t line_between(const vec3_t& p1, const vec3_t& p2) {
    return {p1, p2 - p1};
}


struct plane_t {
    vec3_t pos, normal;
};


struct sphere_t {
    vec3_t pos;
    float r = 0;

    plane_t normal_plane(const vec3_t& loc) const {
        return {loc, (loc - pos).normalized()};
    }
};


struct vec4_t {
    vec3_t a, b;
};

enum struct gshape_type : std::uint32_t {
    plane, sphere
};

struct light_t {
    std::variant<plane_t, sphere_t> obj;
    gshape_type type;
    float strength = 0.0f;
    
    /* function for how quick the light strength should falloff based on distance (x), returns a multiplier */
    std::function<float (float)> falloff = [](float x) -> float {
        /* this is just a step-down function from a to b */
        const float a = 2, b = 10;
        x = std::clamp<float>(x, a, b);
        const float alpha = std::pow(std::numbers::e_v<float>, - (b - a) / (x - a));
        const float beta  = std::pow(std::numbers::e_v<float>, - (b - a) / (b - x));
        return 1 - alpha / (alpha + beta);
    };
};


/* general object type */
enum struct gobj_type : std::uint32_t {
    plane, sphere, light
};

/* general object */
struct gobj_t {
    std::variant<plane_t, sphere_t, light_t> obj;
    gobj_type type;
    float smoothness = 1.0f;
    rgb_t color;
};

vec3_t gobj_get_pos(const gobj_t& gobj, struct notcurses *nc) {
    if (gobj.type == gobj_type::sphere) {
        return std::get<sphere_t>(gobj.obj).pos;
    } else if (gobj.type == gobj_type::plane) {
        return std::get<plane_t>(gobj.obj).pos;
    } else if (gobj.type == gobj_type::light) {
        const auto& light = std::get<light_t>(gobj.obj);
        if (light.type == gshape_type::sphere) {
            return std::get<sphere_t>(light.obj).pos;
        } else if (light.type == gshape_type::plane) {
            return std::get<plane_t>(light.obj).pos;
        }
    }
    notcurses_stop(nc);
    ERR_EXIT("bad gobj.type: %i", gobj.type);
}

void s_intersect(const line_t& line, const sphere_t& sphere, std::vector<float>& out) {
    const float alpha = (line.n.x * line.n.x) + (line.n.y * line.n.y) + (line.n.z * line.n.z);
    const float beta  = 2 * (
        line.n.x * (line.pos.x - sphere.pos.x) +
        line.n.y * (line.pos.y - sphere.pos.y) +
        line.n.z * (line.pos.z - sphere.pos.z));
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
    out.push_back(t1);
    out.push_back(t2);
}

void p_intersect(const line_t& line, const plane_t& plane, std::vector<float>& out) {
    const float alpha = (plane.normal.x * line.n.x + plane.normal.y * line.n.y + plane.normal.z * line.n.z);
    if (alpha == 0) {
        return;
    }
    const float t = (plane.normal.x * (plane.pos.x - line.pos.x) + 
        plane.normal.y * (plane.pos.y - line.pos.y) +
        plane.normal.z * (plane.pos.z - line.pos.z)) / alpha;
    out.push_back(t);
}

/* out is a vector of the times (t) on the line that the intersect occured */
void g_intersect(const line_t& line, const gobj_t& gobj, std::vector<float>& out, struct notcurses *nc) {
    if (gobj.type == gobj_type::sphere) {
        s_intersect(line, std::get<sphere_t>(gobj.obj), out);
    } else if (gobj.type == gobj_type::plane) {
        p_intersect(line, std::get<plane_t>(gobj.obj), out);
    } else if (gobj.type == gobj_type::light) {
        const auto& light = std::get<light_t>(gobj.obj);
        if (light.type == gshape_type::sphere) {
            s_intersect(line, std::get<sphere_t>(light.obj), out);
        } else if (light.type == gshape_type::plane) {
            p_intersect(line, std::get<plane_t>(light.obj), out);
        }
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad gobj.type: %i", gobj.type);
    }
}


vec3_t p_reflect(const vec3_t& vec, const plane_t& plane) {
    return vec - plane.normal * vec.dot(plane.normal) * 2.0f;
}

vec3_t s_reflect(const vec3_t& vec, const sphere_t& sphere, const vec3_t& pos) {
    return p_reflect(vec, sphere.normal_plane(pos));
}

/* reflects vec across the normal plane of gobj at pos */
vec3_t g_reflect(const vec3_t& vec, const gobj_t& gobj, const vec3_t& pos, struct notcurses *nc) {
    if (gobj.type == gobj_type::sphere) {
        return s_reflect(vec, std::get<sphere_t>(gobj.obj), pos);
    } else if (gobj.type == gobj_type::plane) {
        return p_reflect(vec, std::get<plane_t>(gobj.obj));
    } else if (gobj.type == gobj_type::light) {
        const auto& light = std::get<light_t>(gobj.obj);
        if (light.type == gshape_type::sphere) {
            return s_reflect(vec, std::get<sphere_t>(light.obj), pos);
        } else if (light.type == gshape_type::plane) {
            return p_reflect(vec, std::get<plane_t>(light.obj));
        }
    }
    notcurses_stop(nc);
    ERR_EXIT("bad gobj.type: %i", gobj.type);
}


/* using ID = std::uint64_t;

std::uint64_t make_id() {
    static std::random_device device{};
    static std::default_random_engine engine(device());
    static std::uniform_int_distribution<ID> dist;
    return dist(engine);
} */

struct scene_t {
    std::vector<gobj_t*> objects;
    line_t camera_ray;
};

struct char_ex_info_t { /* NOLINT */
    wchar_t ch = L' ';
    rgb_t color{255, 255, 255};

    char_ex_info_t& operator=(const char_ex_info_t& other) = default;

    bool operator==(const char_ex_info_t& other) const {
        return ch == other.ch && color == other.color;
    }

    bool operator!=(const char_ex_info_t& other) const {
        return !(*this == other);
    }
};


std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}



int main() {
    constexpr float fps = 60.0f;
    constexpr std::uint32_t max_light_bounces = 5;

    const std::wstring gradient = L".._,'`^\"-~:;=!><+?|][}{)(\\/trxnuvczijl1fXYUJICLQO0Zmwqpdbkhao*#MW&8%B@$";
    const auto gradient_length = static_cast<std::int32_t>(gradient.length());
    const auto get_gradient = [&gradient, &gradient_length](float x) -> wchar_t {
        return gradient[std::clamp<std::int32_t>(static_cast<std::int32_t>(std::round(x)), 0, gradient_length - 1)];
    };

    const auto wait_us_per_frame = static_cast<std::uint64_t>(1'000'000.0 / fps);


    if (!setlocale(LC_ALL, "")) {
        ERR_EXIT("couldn't set locale");
    }

    /* initialize notcurses */
    notcurses_options opts{};
    opts.flags = NCOPTION_INHIBIT_SETLOCALE;

    struct notcurses *nc = notcurses_init(&opts, nullptr);

    if (nc == nullptr) {
        ERR_EXIT("couldn't initialize notcurses, notcurses_init returned nullptr");
    }

    struct ncplane *plane = notcurses_stdplane(nc);
    ncplane_erase(plane);

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
    const float shift_multiplier = 3.0f;

    /* add objects to scene */
    scene_t scene;
    scene.camera_ray = line_t{vec3_t(0, 0, 0), vec3_t(0, 0, 1)};
    vec3_t& camera_pos = scene.camera_ray.pos;
    vec3_t& camera_n   = scene.camera_ray.n;
    scene.objects.push_back(new gobj_t{
        sphere_t{vec3_t(0, 0, 30), 20.0f},
        gobj_type::sphere,
        1.0f,
        {255, 0, 0}
    });
    /* scene.objects.push_back(new gobj_t{
        plane_t{vec3_t(0, 0, 45), vec3_t(0, 0, -1).normalized()},
        gobj_type::plane,
        0.5f,
        {0, 0, 255}
    }); */
    auto *light_obj = new gobj_t{
        light_t{
            plane_t{vec3_t(-100, 0, 0), vec3_t(1, 0, 0)},
            gshape_type::plane,
            1.0f,
            []([[maybe_unused]] float x) { return 1.0f; }
        },
        gobj_type::light,
        1.0f, /* will not be applied */
        {0, 255, 0}
    };
    auto& light = std::get<light_t>(light_obj->obj);
    /* auto& light_sphere = std::get<sphere_t>(light.obj); */
    scene.objects.push_back(light_obj);
    /* const vec3_t light_orig_pos = light_sphere.pos; */

    std::vector<float> itimes;
    std::uint32_t dimy = 0, dimx = 0;

    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        last_time = get_current_time();

        ncplane_dim_yx(plane, &dimy, &dimx);

        /* revolve the light source */
        /* light_sphere.pos.x = light_orig_pos.x + 15 * std::cos(static_cast<float>(last_time - begin_time) / 1'000'000.0f);
        light_sphere.pos.z = light_orig_pos.z + 15 * std::sin(static_cast<float>(last_time - begin_time) / 1'000'000.0f); */

        
        std::uint32_t chin = notcurses_get_nblock(nc, nullptr);
        if (chin == L'\t') { break; }
        else if (chin == L'q') { camera_pos.y += y_move_speed; }
        else if (chin == L'Q') { camera_pos.y += y_move_speed * shift_multiplier; }
        else if (chin == L'e') { camera_pos.y -= y_move_speed; }
        else if (chin == L'E') { camera_pos.y -= y_move_speed * shift_multiplier; }
        else if (chin == L'w') { camera_pos.z += z_move_speed; }
        else if (chin == L'W') { camera_pos.z += z_move_speed * shift_multiplier; }
        else if (chin == L's') { camera_pos.z -= z_move_speed; }
        else if (chin == L'S') { camera_pos.z -= z_move_speed * shift_multiplier; }
        else if (chin == L'd') { camera_pos.x += x_move_speed; }
        else if (chin == L'D') { camera_pos.x += x_move_speed * shift_multiplier; }
        else if (chin == L'a') { camera_pos.x -= x_move_speed; }
        else if (chin == L'A') { camera_pos.x -= x_move_speed * shift_multiplier; }

        else if (chin == L'j') { camera_n = camera_n.x_rotated(0.1f); }
        else if (chin == L'k') { camera_n = camera_n.x_rotated(-0.1f); }


        /* WARNING: all the following logic basically assumes scene.objects isn't empty */

        for (std::int32_t i = 0; i < dimy; i++) {
            for (std::int32_t j = 0; j < dimx / 2; j++) {
                /* position for i is flipped since notcurses says y down is positive while we want y up is positive */
                vec3_t curpos = vec3_t((j - dimx / 4.0f) * x_mul, (-i + dimy / 2.0f) * y_mul, begin_draw_dist); /* NOLINT */
                line_t ray = line_between(camera_pos, camera_pos + curpos); /* t positive is "forward" */
                char_ex_info_t current_char;

                std::int32_t light_bounces = 0;
                const gobj_t *plight = nullptr;
                const gobj_t *last_obj = nullptr;
                rgb_t original_color;
                line_t line = ray;
                float total_distance = 0.0f;
                float total_smoothness = 1.0f;

                for (; light_bounces < max_light_bounces; light_bounces++) {
                    /* first, determine the closest intersection point out of all objects on the ray */
                    float closest_t = std::numeric_limits<float>::max();
                    const gobj_t *closest_obj = nullptr;
                    vec3_t closest_pos{closest_t, closest_t, closest_t}; /* closest_t is max value right now */
                    /* prgobj = pointer to (potential) reflecting gobj */
                    for (const gobj_t *prgobj : scene.objects) {
                        if (prgobj == last_obj) { continue; }
                        itimes.clear();
                        g_intersect(line, *prgobj, itimes, nc);
                        if (itimes.empty()) { continue; }
                        std::sort(itimes.begin(), itimes.end());
                        float current_t = itimes[0];
                        /* make sure that it is "forward" on the ray, since light has direction */
                        if (current_t < 0.0f) { continue; }
                        float ct = current_t;
                        if (current_t < closest_t) {
                            closest_t = current_t;
                            closest_pos = line.f(current_t);
                            closest_obj = prgobj;
                        }
                    }
                    last_obj = closest_obj;

                    if (closest_obj == nullptr) {
                        break;
                    }
                    if (light_bounces == 0) {
                        original_color = closest_obj->color;
                    }
                    total_distance += (line.pos - closest_pos).mod();
                    if (closest_obj->type == gobj_type::light) {
                        plight = closest_obj;
                        /* done, we are not reflecting off a light */
                        break;
                    }
                    vec3_t newvec = g_reflect(line.n, *closest_obj, closest_pos, nc);
                    total_smoothness *= closest_obj->smoothness;
                    line.n = newvec;
                    line.pos = closest_pos;
                }

                if (plight == nullptr) { /* never got illuminated */
                    if (light_bounces > 0) {
                        current_char = char_ex_info_t{get_gradient(0), original_color};
                    } else {
                        current_char = char_ex_info_t{L' ', {0, 0, 0}};
                    }
                } else {
                    const auto& tlight = std::get<light_t>(plight->obj);
                    float applied_light = tlight.falloff(total_distance) * tlight.strength;
                    rgb_t outcolor;
                    outcolor.x = static_cast<std::int32_t>(static_cast<float>(original_color.x) * applied_light);
                    outcolor.y = static_cast<std::int32_t>(static_cast<float>(original_color.y) * applied_light);
                    outcolor.z = static_cast<std::int32_t>(static_cast<float>(original_color.z) * applied_light);
                    current_char = char_ex_info_t{get_gradient((gradient_length - 1) * applied_light), outcolor};
                }

                ncplane_set_fg_rgb8(plane, current_char.color.x, current_char.color.y, current_char.color.z);
                ncplane_putwc_yx(plane, i, j * 2, current_char.ch);
                ncplane_putwc(plane, current_char.ch);
            }
        }

        ncplane_set_fg_rgb8(plane, 255, 255, 255);
        ncplane_printf_yx(plane, 0, 0, "x: %.2f", camera_pos.x);
        ncplane_printf_yx(plane, 1, 0, "y: %.2f", camera_pos.y);
        ncplane_printf_yx(plane, 2, 0, "z: %.2f", camera_pos.z);

        notcurses_render(nc);
    }

    for (gobj_t *gobj : scene.objects) {
        delete gobj;
    }

    ncplane_destroy(plane);
    notcurses_stop(nc);
    return 0;
}
