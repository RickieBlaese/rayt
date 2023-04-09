#include <algorithm>
#include <chrono>
#include <functional>
#include <limits>
#include <numbers>
#include <random>
#include <variant>
#include <vector>

#include <cinttypes>
#include <cmath>
#include <cstdio>

#include <notcurses/notcurses.h>


#define ERR_EXIT(...) { \
    std::fprintf(stderr, "error: file " __FILE__ ", func %s, line %i: ", __func__, __LINE__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    std::exit(EXIT_FAILURE); \
}


/* general vec3 */
template <typename T>
struct gvec3_t { /* NOLINT */
    T x = 0, y = 0, z = 0;

    gvec3_t() = default;
    gvec3_t(T x, T y, T z) : x(x), y(y), z(z) {}
    gvec3_t(const gvec3_t &other) = default;
    ~gvec3_t() = default;
    bool operator==(const gvec3_t &other) const { return x == other.x && y == other.y && z == other.z; }
    gvec3_t &operator=(const gvec3_t &other) = default;

    gvec3_t operator+(const gvec3_t &other) const { return {x + other.x, y + other.y, z + other.z}; }
    gvec3_t &operator+=(const gvec3_t &other) { return *this = *this + other; }
    gvec3_t operator-() const { return {-x, -y, -z}; }
    gvec3_t operator-(const gvec3_t &other) const { return *this + (-other); }
    gvec3_t &operator-=(const gvec3_t &other) { return *this = *this - other; }

    template <typename Y>
    gvec3_t<T> operator*(const Y &v) const { return {x * v, y * v, z * v}; }
    template <typename Y>
    gvec3_t<T> operator/(const Y &v) const { return {x / v, y / v, z / v}; }

    float mod() const { return std::sqrt(x * x + y * y + z * z); }
    T dot(const gvec3_t &other) const { return x * other.x + y * other.y + z * other.z; }
    gvec3_t normalized() const { return *this / mod(); }
    gvec3_t x_rotated(float theta) const { return {x * std::cos(theta) - z * std::sin(theta), y, x * std::sin(theta) + z * std::cos(theta)}; }
    gvec3_t y_rotated(float theta) const { return {x, y * std::cos(theta) - z * std::sin(theta), y * std::sin(theta) + z * std::cos(theta)}; }
};

using vec3_t = gvec3_t<float>;
using rgb_t = gvec3_t<std::int32_t>;


/* markiplier -\ 
 *             |
 *             v */
rgb_t multiplier(const rgb_t &original_color, float k) {
    rgb_t outcolor;
    outcolor.x = static_cast<std::int32_t>(static_cast<float>(original_color.x) * k);
    outcolor.y = static_cast<std::int32_t>(static_cast<float>(original_color.y) * k);
    outcolor.z = static_cast<std::int32_t>(static_cast<float>(original_color.z) * k);
    return outcolor;
}

rgb_t average_colors(const rgb_t &a, const rgb_t &b) {
    return multiplier(a + b,  0.5f);
}

inline void sort_by_dist(std::vector<vec3_t> &vecs, const vec3_t &pos) {
    std::sort(vecs.begin(), vecs.end(), [&pos](const vec3_t &a, const vec3_t &b) { return (pos - a).mod() < (pos - b).mod(); });
}


struct line_t {
    vec3_t pos, n;

    vec3_t f(float t) const {
        return {pos.x + n.x * t, pos.y + n.y * t, pos.z + n.z * t};
    }
};

line_t line_between(const vec3_t &p1, const vec3_t &p2) {
    return {p1, p2 - p1};
}


struct plane_t {
    vec3_t pos, normal;
};


struct bounded_plane_t {
    plane_t plane; 
    /* a, b is begin coords with respect to plane.pos
     * c, d are "length" or "width" of a, b respectively
     * c, d must be positive
     * if normal axis is x then these are y, z coords
     * if normal axis is y then these are x, z coords
     * if normal axis is z then these are x, y coords */
    float a = 0, b = 0, c = 0, d = 0;
};

struct sphere_t {
    vec3_t pos;
    float r = 0;

    plane_t normal_plane(const vec3_t &loc) const {
        return {loc, (loc - pos).normalized()};
    }
};

/* general variant indices, to be compared with variant.index() */
constexpr std::size_t gtype_sphere = 0;
constexpr std::size_t gtype_plane = 1;
constexpr std::size_t gtype_bounded_plane = 2;
constexpr std::size_t gtype_light = 3;

struct light_t {
    std::variant<sphere_t, plane_t, bounded_plane_t> obj;
    
    /* function for how quick the light strength should falloff based on distance (x) */
    std::function<float (float)> strength = [](float x) -> float {
        /* this is just a step-down function from a to b */
        const float a = 0, b = 200;
        x = std::clamp<float>(x, a, b);
        const float alpha = std::pow(std::numbers::e_v<float>, - (b - a) / (x - a));
        const float beta  = std::pow(std::numbers::e_v<float>, - (b - a) / (b - x));
        return 1.0f - alpha / (alpha + beta);
    };
};


/* general object */
struct gobj_t {
    std::variant<sphere_t, plane_t, bounded_plane_t, light_t> obj;
    rgb_t color;
    bool mirror = false;
};

enum struct axis_t : std::uint32_t {
    x, y, z
};

axis_t normal_to_axis(const vec3_t &normal, struct notcurses *nc) {
    if (normal.x == 0 && normal.z == 0) {
        return axis_t::y;
    } else if (normal.x == 0 && normal.y == 0) {
        return axis_t::z;
    } else if (normal.y == 0 && normal.z == 0) {
        return axis_t::x;
    }
    notcurses_stop(nc);
    ERR_EXIT("bad nontrivial normal passed: vec3_t(%.2f, %.2f, %.2f)", normal.x, normal.y, normal.z);
}

vec3_t gobj_get_pos(const gobj_t &gobj, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        return std::get<sphere_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_plane) {
        return std::get<plane_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_bounded_plane) {
        return std::get<bounded_plane_t>(gobj.obj).plane.pos;
    } else if (gobj.obj.index() == gtype_light) {
        const auto &light = std::get<light_t>(gobj.obj);
        if (light.obj.index() == gtype_sphere) {
            return std::get<sphere_t>(light.obj).pos;
        } else if (light.obj.index() == gtype_plane) {
            return std::get<plane_t>(light.obj).pos;
        } else if (light.obj.index() == gtype_bounded_plane) {
            return std::get<bounded_plane_t>(light.obj).plane.pos;
        }
        notcurses_stop(nc);
        ERR_EXIT("bad light variant index: std::get<light_t>(gobj.obj).obj.index() = %zu", light.obj.index());
    }
    notcurses_stop(nc);
    ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
}


void s_intersect(const line_t &line, const sphere_t &sphere, std::vector<float> &out) {
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

void p_intersect(const line_t &line, const plane_t &plane, std::vector<float> &out) {
    const float alpha = (plane.normal.x * line.n.x + plane.normal.y * line.n.y + plane.normal.z * line.n.z);
    if (alpha == 0) {
        return;
    }
    const float t = (plane.normal.x * (plane.pos.x - line.pos.x) + 
        plane.normal.y * (plane.pos.y - line.pos.y) +
        plane.normal.z * (plane.pos.z - line.pos.z)) / alpha;
    out.push_back(t);
}

void bp_intersect(const line_t &line, const bounded_plane_t &bounded_plane, std::vector<float> &out, struct notcurses *nc) {
    static std::vector<float> thisout;
    thisout.clear();
    p_intersect(line, bounded_plane.plane, thisout);
    axis_t axis = normal_to_axis(bounded_plane.plane.normal, nc);
    switch (axis) {
        case axis_t::x:
            for (float t : thisout) {
                const vec3_t pos = line.f(t);
                /* is bad / outside */
                if ((pos.y - bounded_plane.plane.pos.y < bounded_plane.a || pos.y - bounded_plane.plane.pos.y > bounded_plane.a + bounded_plane.c) ||
                    (pos.z - bounded_plane.plane.pos.z < bounded_plane.b || pos.z - bounded_plane.plane.pos.z > bounded_plane.b + bounded_plane.d)) {
                    continue;
                }
                out.push_back(t);
            }
            break;
        case axis_t::y:
            for (float t : thisout) {
                const vec3_t pos = line.f(t);
                /* is bad / outside */
                if ((pos.x - bounded_plane.plane.pos.x < bounded_plane.a || pos.x - bounded_plane.plane.pos.x > bounded_plane.a + bounded_plane.c) ||
                    (pos.z - bounded_plane.plane.pos.z < bounded_plane.b || pos.z - bounded_plane.plane.pos.z > bounded_plane.b + bounded_plane.d)) {
                    continue;
                }
                out.push_back(t);
            }
            break;
        case axis_t::z:
            for (float t : thisout) {
                const vec3_t pos = line.f(t);
                /* is bad / outside */
                if ((pos.x - bounded_plane.plane.pos.x < bounded_plane.a || pos.x - bounded_plane.plane.pos.x > bounded_plane.a + bounded_plane.c) ||
                    (pos.y - bounded_plane.plane.pos.y < bounded_plane.b || pos.y - bounded_plane.plane.pos.y > bounded_plane.b + bounded_plane.d)) {
                    continue;
                }
                out.push_back(t);
            }
            break;
    }
}

/* out is a vector of the times (t) on the line that the intersect occured */
void g_intersect(const line_t &line, const gobj_t &gobj, std::vector<float> &out, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        s_intersect(line, std::get<sphere_t>(gobj.obj), out);
    } else if (gobj.obj.index() == gtype_plane) {
        p_intersect(line, std::get<plane_t>(gobj.obj), out);
    } else if (gobj.obj.index() == gtype_bounded_plane) {
        bp_intersect(line, std::get<bounded_plane_t>(gobj.obj), out, nc);
    } else if (gobj.obj.index() == gtype_light) {
        const auto &light = std::get<light_t>(gobj.obj);
        if (light.obj.index() == gtype_sphere) {
            s_intersect(line, std::get<sphere_t>(light.obj), out);
        } else if (light.obj.index() == gtype_plane) {
            p_intersect(line, std::get<plane_t>(light.obj), out);
        } else if (light.obj.index() == gtype_bounded_plane) {
            bp_intersect(line, std::get<bounded_plane_t>(light.obj), out, nc);
        } else {
            notcurses_stop(nc);
            ERR_EXIT("bad light variant index: std::get<light_t>(gobj.obj).obj.index() = %zu", light.obj.index());
        }
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
}


vec3_t p_reflect(const vec3_t &vec, const plane_t &plane) {
    return vec - plane.normal.normalized() * vec.dot(plane.normal.normalized()) * 2.0f;
}


vec3_t s_reflect(const vec3_t &vec, const sphere_t &sphere, const vec3_t &pos) {
    return p_reflect(vec, sphere.normal_plane(pos));
}

/* reflects vec across the normal plane of gobj at pos */
vec3_t g_reflect(const vec3_t &vec, const gobj_t &gobj, const vec3_t &pos, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        return s_reflect(vec, std::get<sphere_t>(gobj.obj), pos);
    } else if (gobj.obj.index() == gtype_plane) {
        return p_reflect(vec, std::get<plane_t>(gobj.obj));
    } else if (gobj.obj.index() == gtype_bounded_plane) {
        return p_reflect(vec, std::get<bounded_plane_t>(gobj.obj).plane);
    } else if (gobj.obj.index() == gtype_light) {
        const auto &light = std::get<light_t>(gobj.obj);
        if (light.obj.index() == gtype_sphere) {
            return s_reflect(vec, std::get<sphere_t>(light.obj), pos);
        } else if (light.obj.index() == gtype_plane) {
            return p_reflect(vec, std::get<plane_t>(light.obj));
        } else if (light.obj.index() == gtype_bounded_plane) {
            return p_reflect(vec, std::get<bounded_plane_t>(light.obj).plane);
        }
        notcurses_stop(nc);
        ERR_EXIT("bad light variant index: std::get<light_t>(gobj.obj).obj.index() = %zu", light.obj.index());
    }
    notcurses_stop(nc);
    ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
}


/* rotate vec around pos */
vec3_t rotated_x_about(const vec3_t &vec, const vec3_t &pos, float theta) {
    return (vec - pos).x_rotated(theta) + pos;
}

/* rotate vec around pos */
vec3_t rotated_y_about(const vec3_t &vec, const vec3_t &pos, float theta) {
    return (vec - pos).y_rotated(theta) + pos;
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

    char_ex_info_t &operator=(const char_ex_info_t &other) = default;

    bool operator==(const char_ex_info_t &other) const {
        return ch == other.ch && color == other.color;
    }

    bool operator!=(const char_ex_info_t &other) const {
        return !(*this == other);
    }
};


std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}



int main() {
    constexpr float fps = 60.0f;
    constexpr std::uint32_t max_light_bounces = 20;

    const std::wstring gradient = L" ._,'`^\"-~:;=!><+?|][}{)(\\/trxnuovczmwaihqpdbkfjl1XYFGHNUJICLQO0Z#MW&8%B@$";
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
    const float begin_draw_dist = 50.0f;

    /* aspect ratio ? */
    const float x_mul = 1.0f;
    const float y_mul = 1.0f;

    const float z_move_speed = 1.0f;
    const float x_move_speed = 1.0f;
    const float y_move_speed = 1.0f;
    const float shift_multiplier = 3.0f;
    const float angle_increment = std::numbers::pi / 30.0f;

    const float minimum_color_multiplier = 0.1f;


    /* --- add objects to scene --- */
    scene_t scene;
    scene.camera_ray = line_t{vec3_t(0, 0, 0), vec3_t(0, 0, begin_draw_dist)};
    vec3_t &camera_pos = scene.camera_ray.pos;
    vec3_t &camera_n   = scene.camera_ray.n;

    
    /* spheres */
    scene.objects.push_back(new gobj_t{
        sphere_t{vec3_t(0, 0, 30), 20.0f},
        {0, 255, 0}
    });
    scene.objects.push_back(new gobj_t{
        sphere_t{vec3_t(0, 0, 80), 20.0f},
        {0, 0, 255}
    });


    /* planar mirrors */
    scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(0, 0, 105), vec3_t(0, 0, -1)},
            -200, -200, 400, 400
        },
        {0, 0, 0},
        true,
    });
    /* scene.objects.push_back(new gobj_t{
        plane_t{vec3_t(0, 0, 5), vec3_t(0, 0, 1)},
        {0, 0, 0},
        true,
    }); */


    /* bounded planes */
    scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(0, 0, 1), vec3_t(0, 0, -1)},
            -4, -4, 8, 8
        },
        {255, 255, 255}
    });


    /* lights */
    scene.objects.push_back(new gobj_t{
        light_t{
            bounded_plane_t{
                plane_t{vec3_t(100, 0, 0), vec3_t(-1, 0, 0)},
                -200, -200, 400, 400
            },
            []([[maybe_unused]] float x) { return std::pow(std::numbers::e_v<float>, -x/100.0f) + 0.3f; }
        },
        {255, 127, 127}
    });
    scene.objects.push_back(new gobj_t{
        light_t{
            bounded_plane_t{
                plane_t{vec3_t(-100, 0, 0), vec3_t(1, 0, 0)},
                -200, -200, 400, 400
            },
            []([[maybe_unused]] float x) { return std::pow(std::numbers::e_v<float>, -x/100.0f) + 0.3f; }
        },
        {255, 127, 127}
    });
    scene.objects.push_back(new gobj_t{
        light_t{
            bounded_plane_t{
                plane_t{vec3_t(0, 0, -30), vec3_t(0, 0, 1)},
                -200, -200, 400, 400
            },
            []([[maybe_unused]] float x) { return std::pow(std::numbers::e_v<float>, -x/100.0f) + 0.3f; }
        },
        {255, 127, 127}
    });
    scene.objects.push_back(new gobj_t{
        light_t{
            bounded_plane_t{
                plane_t{vec3_t(0, 100, 0), vec3_t(0, -1, 0)},
                -200, -200, 400, 400
            },
            []([[maybe_unused]] float x) { return std::pow(std::numbers::e_v<float>, -x/100.0f) + 0.3f; }
        },
        {255, 127, 127}
    });
    scene.objects.push_back(new gobj_t{
        light_t{
            bounded_plane_t{
                plane_t{vec3_t(0, -100, 0), vec3_t(0, 1, 0)},
                -200, -200, 400, 400
            },
            []([[maybe_unused]] float x) { return std::pow(std::numbers::e_v<float>, -x/100.0f) + 0.3f; }
        },
        {255, 127, 127}
    });

    std::vector<float> itimes;
    std::uint32_t dimy = 0, dimx = 0;
    float camera_rotation_x = 0.0f; /* accounts for both x and z */
    float camera_rotation_y = 0.0f;

    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        last_time = get_current_time();

        ncplane_dim_yx(plane, &dimy, &dimx);
        
        std::uint32_t chin = notcurses_get_nblock(nc, nullptr);
        vec3_t move_d;
        if (chin == L'\t') { break; }
        else if (chin == L'q') { move_d = {0, y_move_speed, 0}; }
        else if (chin == L'Q') { move_d = {0, y_move_speed * shift_multiplier, 0}; }
        else if (chin == L'e') { move_d = {0, -y_move_speed, 0}; }
        else if (chin == L'E') { move_d = {0, -y_move_speed * shift_multiplier, 0}; }
        else if (chin == L'w') { move_d = {0, 0, z_move_speed}; }
        else if (chin == L'W') { move_d = {0, 0, z_move_speed * shift_multiplier}; }
        else if (chin == L's') { move_d = {0, 0, -z_move_speed}; }
        else if (chin == L'S') { move_d = {0, 0, -z_move_speed * shift_multiplier}; }
        else if (chin == L'd') { move_d = {x_move_speed, 0, 0}; }
        else if (chin == L'D') { move_d = {x_move_speed * shift_multiplier, 0, 0}; }
        else if (chin == L'a') { move_d = {-x_move_speed, 0, 0}; }
        else if (chin == L'A') { move_d = {-x_move_speed * shift_multiplier, 0, 0}; }

        else if (chin == L'j') { camera_n = rotated_x_about(camera_n, camera_pos,  angle_increment); camera_rotation_x += angle_increment; }
        else if (chin == L'l') { camera_n = rotated_x_about(camera_n, camera_pos, -angle_increment); camera_rotation_x -= angle_increment; }
        else if (chin == L'i' && camera_rotation_y < std::numbers::pi / 2.0f - angle_increment) {
            camera_n = rotated_y_about(camera_n, camera_pos, angle_increment);
            camera_rotation_y += angle_increment;
        } else if (chin == L'k' && camera_rotation_y > - (std::numbers::pi / 2.0f - angle_increment)) {
            camera_n = rotated_y_about(camera_n, camera_pos, -angle_increment);
            camera_rotation_y -= angle_increment;
        }

        move_d = move_d.x_rotated(camera_rotation_x);
        camera_pos += move_d;
        camera_n   += move_d;


        /* WARNING: all the following logic basically assumes scene.objects isn't empty */

        for (std::int32_t i = 0; i < dimy; i++) {
            for (std::int32_t j = 0; j < dimx / 2; j++) {
                /* position for i is flipped since notcurses says y down is positive while we want y up is positive */
                vec3_t curpos = vec3_t((j - dimx / 4.0f) * x_mul, (-i + dimy / 2.0f) * y_mul, 0); /* NOLINT */
                curpos = curpos.y_rotated(camera_rotation_y).x_rotated(camera_rotation_x); /* is still "at the origin" */
                line_t ray = line_between(camera_pos, camera_n + curpos); /* t positive is "forward" */
                char_ex_info_t current_char;

                std::int32_t light_bounces = 0;
                const gobj_t *original_obj = nullptr;
                const gobj_t *light = nullptr;
                const gobj_t *last_obj = nullptr;
                line_t line = ray;
                float total_distance = 0.0f;

                for (light_bounces = 0; light_bounces < max_light_bounces; light_bounces++) {
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
                        if (current_t < closest_t) {
                            closest_t = current_t;
                            closest_pos = line.f(current_t);
                            closest_obj = prgobj;
                        }
                    }

                    if (closest_obj == nullptr) {
                        break;
                    }
                    last_obj = closest_obj;

                    if (original_obj == nullptr && !closest_obj->mirror) {
                        original_obj = closest_obj;
                    }

                    total_distance += (line.pos - closest_pos).mod();
                    if (closest_obj->obj.index() == gtype_light) {
                        /* done, we are not reflecting off a light */
                        light = closest_obj;
                        break;
                    }
                    vec3_t newvec = g_reflect(line.n, *closest_obj, closest_pos, nc);
                    line.n = newvec;
                    line.pos = closest_pos;
                }


                rgb_t color;
                if (original_obj != nullptr) {
                    color = original_obj->color; /* * !original_obj->mirror + last_obj->color * original_obj->mirror *//* funny branchless */
                    if (original_obj->mirror) {
                        color = last_obj->color;
                    }
                }

                current_char = char_ex_info_t{L' ', {0, 0, 0}};
                if (last_obj != nullptr) {
                    if (light == nullptr) { /* never got illuminated */
                        if (light_bounces > 0) {
                            current_char = char_ex_info_t{L'.', multiplier(color, minimum_color_multiplier)};
                        }
                    } else {
                        const auto &tlight = std::get<light_t>(light->obj);
                        float applied_light = tlight.strength(total_distance);
                        if (light_bounces > 0) {
                            color.x = static_cast<std::int32_t>(static_cast<float>(color.x * light->color.x) / 255.0f);
                            color.y = static_cast<std::int32_t>(static_cast<float>(color.y * light->color.y) / 255.0f);
                            color.z = static_cast<std::int32_t>(static_cast<float>(color.z * light->color.z) / 255.0f);
                        }
                        current_char = char_ex_info_t{get_gradient(static_cast<float>(gradient_length - 1) * applied_light), multiplier(color, applied_light)};
                    }
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
        ncplane_printf_yx(plane, 3, 0, "nx: %.2f", camera_n.x);
        ncplane_printf_yx(plane, 4, 0, "ny: %.2f", camera_n.y);
        ncplane_printf_yx(plane, 5, 0, "nz: %.2f", camera_n.z);

        notcurses_render(nc);
    }

    for (gobj_t *gobj : scene.objects) {
        delete gobj;
    }

    ncplane_destroy(plane);
    notcurses_stop(nc);
    return 0;
}
