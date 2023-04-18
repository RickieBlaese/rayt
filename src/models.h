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
    gvec3_t<T> operator*(const Y &v) const { return {x * v, y * v, z * v}; }

    template <typename Y>
    gvec3_t<T> operator/(const Y &v) const { return {x / v, y / v, z / v}; }

    float mod() const { return std::sqrt(x * x + y * y + z * z); }
    T dot(const gvec3_t<T> &other) const { return x * other.x + y * other.y + z * other.z; }
    gvec3_t<T> normalized() const { return *this / mod(); }
    gvec3_t<T> x_rotated(float theta) const { return {x * std::cos(theta) - z * std::sin(theta), y, x * std::sin(theta) + z * std::cos(theta)}; }
    gvec3_t<T> y_rotated(float theta) const { return {x, y * std::cos(theta) - z * std::sin(theta), y * std::sin(theta) + z * std::cos(theta)}; }
};

using vec3_t = gvec3_t<float>;
using rgb_t = gvec3_t<std::int32_t>;

/* markiplier -\ 
 *             |
 *             v */
rgb_t multiplier(const rgb_t &original_color, float k);

rgb_t average_colors(const rgb_t &a, const rgb_t &b);

inline void sort_by_dist(std::vector<vec3_t> &vecs, const vec3_t &pos);


struct line_t {
    vec3_t pos, n;

    vec3_t f(float t) const;
};

line_t line_between(const vec3_t &p1, const vec3_t &p2);


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

    plane_t normal_plane(const vec3_t &loc) const;
};

struct rect_t {
    vec3_t pos;
    bounded_plane_t bottom, top, left, right, back, front;
};

/* size.x, size.y, size.z are width in that direction, should be positive 
 * pos is the origin position, this should be the minimum
 * coords of the cube, i.e. bottom left back vertex */
rect_t create_rect(bool normal_outward, const vec3_t &pos, const vec3_t &size);


/* general variant indices, to be compared with variant.index() */
constexpr std::size_t gtype_sphere = 0;
constexpr std::size_t gtype_plane = 1;
constexpr std::size_t gtype_bounded_plane = 2;

using shape_variant_t = std::variant<sphere_t, plane_t, bounded_plane_t>;

float default_light_strength(float x);

/* general object */
struct gobj_t {
    shape_variant_t obj;
    rgb_t color;
    bool mirror = false;
    bool light = false;
    
    /* function for how quick the light strength should falloff based on distance (x) */
    float (*strength)(float) = default_light_strength;
};

enum struct axis_t : std::uint32_t {
    x, y, z
};

axis_t normal_to_axis(const vec3_t &normal, struct notcurses *nc);

vec3_t gobj_get_pos(const gobj_t &gobj, struct notcurses *nc);


void s_intersect(const line_t &line, const sphere_t &sphere, std::vector<float> &out);

void p_intersect(const line_t &line, const plane_t &plane, std::vector<float> &out);

void bp_intersect(const line_t &line, const bounded_plane_t &bounded_plane, std::vector<float> &out, struct notcurses *nc);

/* out is a vector of the times (t) on the line that the intersect occured */
void g_intersect(const line_t &line, const gobj_t &gobj, std::vector<float> &out, struct notcurses *nc);

vec3_t p_reflect(const vec3_t &vec, const plane_t &plane);

vec3_t s_reflect(const vec3_t &vec, const sphere_t &sphere, const vec3_t &pos);

/* reflects vec across the nomal plane of gobj at pos */
vec3_t g_reflect(const vec3_t &vec, const gobj_t &gobj, const vec3_t &pos, struct notcurses *nc);

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
    float minimum_color_multiplier = 0.0f;
    std::wstring gradient;

    wchar_t get_gradient(float x);
    void intersect_ray(const line_t &line, const gobj_t *last_obj, gobj_t *&closest_obj, float &closest_t, struct notcurses *nc);
    void render_ray(const line_t &line, rgb_t &outcolor, wchar_t &outchar, struct notcurses *nc);
};

/* allocates gobj_t */
void add_rect_light(scene_t &scene, const rect_t &rect, const rgb_t &color, bool mirror, const decltype(gobj_t::strength)& strength);

void add_rect(scene_t &scene, const rect_t &rect, const rgb_t &color, bool mirror);


#endif
