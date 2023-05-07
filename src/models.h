#ifndef MODELS_H
#define MODELS_H

#include "common.h"

/* general vec3 */
template <typename T>
struct gvec3_t { /* NOLINT */
    T x = 0, y = 0, z = 0;

    gvec3_t() = default;
    gvec3_t(T x, T y, T z) : x(x), y(y), z(z) {}
    gvec3_t(const gvec3_t<T> &other) = default;
    ~gvec3_t() = default;
    bool operator==(const gvec3_t<T> &other) const { return x == other.x && y == other.y && z == other.z; }
    gvec3_t<T> &operator=(const gvec3_t<T> &other) = default;

    gvec3_t<T> operator+(const gvec3_t<T> &other) const { return {x + other.x, y + other.y, z + other.z}; }
    gvec3_t<T> &operator+=(const gvec3_t<T> &other) { return *this = *this + other; }
    gvec3_t<T> operator-() const { return {-x, -y, -z}; }
    gvec3_t<T> operator-(const gvec3_t<T> &other) const { return *this + (-other); }
    gvec3_t<T> &operator-=(const gvec3_t<T> &other) { return *this = *this - other; }

    template <typename Y>
    [[nodiscard]] gvec3_t<T> operator*(const Y &v) const { return {x * v, y * v, z * v}; }

    template <typename Y>
    [[nodiscard]] gvec3_t<T> operator/(const Y &v) const { return {x / v, y / v, z / v}; }

    [[nodiscard]] float mod() const { return std::sqrt(x * x + y * y + z * z); }
    [[nodiscard]] T dot(const gvec3_t<T> &other) const { return x * other.x + y * other.y + z * other.z; }
    [[nodiscard]] gvec3_t<T> cross(const gvec3_t<T> &other) const { return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x}; }
    [[nodiscard]] gvec3_t<T> normalized() const { return *this / mod(); }
    [[nodiscard]] gvec3_t<T> x_rotated(float theta) const { return {x * std::cos(theta) - z * std::sin(theta), y, x * std::sin(theta) + z * std::cos(theta)}; }
    [[nodiscard]] gvec3_t<T> y_rotated(float theta) const { return {x, y * std::cos(theta) - z * std::sin(theta), y * std::sin(theta) + z * std::cos(theta)}; }
};

using vec3_t = gvec3_t<float>;
using rgb_t = gvec3_t<std::int32_t>;

/* markiplier -\ 
 *             |
 *             v */
rgb_t multiplier(const rgb_t &original_color, float k);

rgb_t color_multiplier(const rgb_t &a, const rgb_t &b);

rgb_t average_colors(const rgb_t &a, const rgb_t &b);

rgb_t clamp(const rgb_t &color);

inline void sort_by_dist(std::vector<vec3_t> &vecs, const vec3_t &pos);


struct line_t {
    vec3_t pos, n;

    vec3_t f(float t) const;
};

line_t line_between(const vec3_t &p1, const vec3_t &p2);


struct plane_t {
    vec3_t pos, normal;
};


/* convex only
 * assumes they already lie on a single plane */
template <std::uint32_t C> requires std::ratio_greater_equal<std::ratio<C>, std::ratio<3>>::value
struct polygon_t {
    vec3_t v[C];

    vec3_t get_pos() const {
        vec3_t acc;
        for (std::uint32_t i = 0; i < C; i++) {
            acc += v[i];
        }
        return acc / static_cast<float>(C);
    }

    plane_t plane() const {
        return {get_pos(), (v[1] - v[0]).cross(v[2] - v[1]).normalized()};
    }

    /* reversed winding order, i.e. flipped orientation */
    polygon_t<C> reversed_wo() const {
        polygon_t<C> r;
        std::reverse_copy(std::begin(v), std::end(v), std::begin(r.v));
        return r;
    }
};


struct sphere_t {
    vec3_t pos;
    float r = 0;

    plane_t normal_plane(const vec3_t &loc) const;
};

using rect_t = polygon_t<4>;

using triangle_t = polygon_t<3>;

struct rectprism_t {
    rect_t bottom, top, left, right, back, front;

    rectprism_t reversed_wo() const;
    vec3_t get_pos() const;
};

/* size.x, size.y, size.z are width in that directprismion, should be positive 
 * pos is the origin position, this should be the minimum
 * coords of the cube, i.e. bottom left back vertex */
rectprism_t create_rectprism(const vec3_t &pos, const vec3_t &size);

struct cylinder_t {
    line_t l;
    float r = 0;

    plane_t normal_plane(const vec3_t &loc) const;
};

/* general variant indices, to be compared with variant.index() */
constexpr std::size_t gtype_sphere = 0;
constexpr std::size_t gtype_plane = 1;
constexpr std::size_t gtype_rect = 2;
constexpr std::size_t gtype_triangle = 3;
constexpr std::size_t gtype_cylinder = 4;

using shape_variant_t = std::variant<sphere_t, plane_t, rect_t, triangle_t, cylinder_t>;

float default_light_strength([[maybe_unused]] float x);

struct gobj_t;

rgb_t default_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos);
float default_opacity_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos);
float default_roughness_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos);

struct scene_t;

/* general object */
struct gobj_t {
    shape_variant_t obj;
    rgb_t color;
    bool mirror = false;
    bool light = false;
    
    /* function for how quick the light strength should falloff based on distance (x) */
    float (*strength)([[maybe_unused]] float x) = default_light_strength;
    /* returns color based on world position */
    rgb_t (*texture)([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) = default_texture;
    float (*opacity_texture)([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) = default_opacity_texture;
    float (*roughness_texture)([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) = default_roughness_texture;

    scene_t *portal = nullptr;
    float roughness = 0.0f;
    bool transparent = false;
    float opacity = 1.0f;
    bool hidden = false;
};

vec3_t gobj_get_pos(const gobj_t &gobj, struct notcurses *nc);

void gobj_set_pos(gobj_t &gobj, const vec3_t &pos, struct notcurses *nc);


bool s_intersect(const line_t &line, const sphere_t &sphere, std::pair<std::optional<float>, std::optional<float>> &ptimes);

bool p_intersect(const line_t &line, const plane_t &plane, std::pair<std::optional<float>, std::optional<float>> &ptimes);

template <std::uint32_t C>
bool pg_intersect(const line_t &line, const polygon_t<C> &polygon, std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc);

bool cl_intersect(const line_t &line, const cylinder_t &cylinder, std::pair<std::optional<float>, std::optional<float>> &ptimes);

bool g_intersect(const line_t &line, const gobj_t &gobj, std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc);

vec3_t p_reflect(const vec3_t &vec, const plane_t &plane, vec3_t &normal);

/* reflects vec across the nomal plane of gobj at pos */
vec3_t g_reflect(const vec3_t &vec, const gobj_t &gobj, const vec3_t &pos, struct notcurses *nc, vec3_t &normal);

/* rotate vec around pos */
vec3_t rotated_x_about(const vec3_t &vec, const vec3_t &pos, float theta);

/* rotate vec around pos */
vec3_t rotated_y_about(const vec3_t &vec, const vec3_t &pos, float theta);


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
    std::uint32_t max_light_bounces = 20;
    std::uint32_t samples_per_ray = 10;
    float minimum_color_multiplier = 0.01f;
    std::wstring gradient;

    wchar_t get_gradient(float x);
    void intersect_ray(const line_t &line, const gobj_t *last_obj, gobj_t *&closest_obj, float &closest_t, struct notcurses *nc);
    /* returns how many times it bounced light */
    std::uint64_t render_ray(const line_t &ray, rgb_t &outcolor, wchar_t &outchar, float &applied_light, struct notcurses *nc);
};

struct world_t {
    std::vector<scene_t*> scenes;
};

/* allocates gobj_t */
void add_rectprism_light(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror, const decltype(gobj_t::strength) &strength);

/* allocates gobj_t */
void add_rectprism(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror, float roughness);


struct intersection_t {
    rgb_t color;
    float roughness = 0.0f, dist = 0.0f,
          c = 0.0f, p = 1.0f, brdf = 0.0f;
    bool light = false;
    bool transparent = false;
    float opacity = 0.0f;
};

#endif
