#include "models.h"


rgb_t multiplier(const rgb_t &original_color, float k) {
    rgb_t outcolor;
    outcolor.x = static_cast<std::int32_t>(static_cast<float>(original_color.x) * k);
    outcolor.y = static_cast<std::int32_t>(static_cast<float>(original_color.y) * k);
    outcolor.z = static_cast<std::int32_t>(static_cast<float>(original_color.z) * k);
    return outcolor;
}

rgb_t average_colors(const rgb_t &a, const rgb_t &b) {
    return multiplier(a + b, 0.5f);
}

inline void sort_by_dist(std::vector<vec3_t> &vecs, const vec3_t &pos) {
    std::sort(vecs.begin(), vecs.end(), [&pos](const vec3_t &a, const vec3_t &b) { return (pos - a).mod() < (pos - b).mod(); });
}

vec3_t line_t::f(float t) const {
    return {pos.x + n.x * t, pos.y + n.y * t, pos.z + n.z * t};
}

line_t line_between(const vec3_t &p1, const vec3_t &p2) {
    return {p1, p2 - p1};
}

plane_t sphere_t::normal_plane(const vec3_t &loc) const {
    /* multiplying by r accounts for inverse spheres, making the normal opposite if r is negative */
    return {loc, ((loc - pos)).normalized()};
}

rectprism_t create_rectprism(const vec3_t &pos, const vec3_t &size) {
    return {
        vec3_t(pos.x + size.x / 2.0f, pos.y + size.y / 2.0f, pos.z + size.z / 2.0f),
        /* bottom */
        rect_t{{
            pos, pos + vec3_t(0, 0, size.z),
            pos + vec3_t(size.x, 0, size.z), pos + vec3_t(size.x, 0, 0)
        }},
        /* top */
        rect_t{{
            pos + vec3_t(0, size.y, 0), pos + vec3_t(0, size.y, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, size.y, 0)
        }},
        /* left */
        rect_t{{
            pos, pos + vec3_t(0, 0, size.z),
            pos + vec3_t(0, size.y, size.z), pos + vec3_t(0, size.y, 0)
        }},
        /* right */
        rect_t{{
            pos + vec3_t(size.x, 0, 0), pos + vec3_t(size.x, 0, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, size.y, 0)
        }},
        /* back */
        rect_t{{
            pos, pos + vec3_t(0, size.y, 0),
            pos + vec3_t(size.x, size.y, 0), pos + vec3_t(size.x, 0, 0)
        }},
        /* front */
        rect_t{{
            pos, pos + vec3_t(0, size.y, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, 0, size.z)
        }}
    };
}

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


float default_light_strength([[maybe_unused]] float x) {
    /* this is just a smooth step-down function from a to b */
    const float a = 0, b = 200;
    x = std::clamp<float>(x, a, b);
    const float alpha = std::pow(std::numbers::e_v<float>, - (b - a) / (x - a));
    const float beta  = std::pow(std::numbers::e_v<float>, - (b - a) / (b - x));
    return 1.0f - alpha / (alpha + beta);
}

rgb_t default_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.color;
}

float default_opacity_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.opacity;
}

float default_roughness_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.roughness;
}

vec3_t gobj_get_pos(const gobj_t &gobj, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        return std::get<sphere_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_plane) {
        return std::get<plane_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_rect) {
        return std::get<rect_t>(gobj.obj).plane().pos;
    } else if (gobj.obj.index() == gtype_triangle) {
        return std::get<triangle_t>(gobj.obj).plane().pos;
    }
    notcurses_stop(nc);
    ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
}

void gobj_set_pos(gobj_t &gobj, const vec3_t &pos, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        std::get<sphere_t>(gobj.obj).pos = pos;
    } else if (gobj.obj.index() == gtype_plane) {
        std::get<plane_t>(gobj.obj).pos = pos;
    } else if (gobj.obj.index() == gtype_rect) {
        auto &rect = std::get<rect_t>(gobj.obj);
        const vec3_t d = pos - rect.plane().pos;
        rect.v[0] = rect.v[0] + d;
        rect.v[1] = rect.v[1] + d;
        rect.v[2] = rect.v[2] + d;
        rect.v[3] = rect.v[3] + d;
    } else if (gobj.obj.index() == gtype_triangle) {
        auto &triangle = std::get<triangle_t>(gobj.obj);
        const vec3_t d = pos - triangle.plane().pos;
        triangle.v[0] = triangle.v[0] + d;
        triangle.v[1] = triangle.v[1] + d;
        triangle.v[2] = triangle.v[2] + d;
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
}

bool s_intersect(const line_t &line, const sphere_t &sphere, std::pair<std::optional<float>, std::optional<float>> &ptimes) {
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
        return false;
    }
    ptimes.first  = (-beta + std::sqrt(discr)) / (2 * alpha);
    ptimes.second = (-beta - std::sqrt(discr)) / (2 * alpha);
    return true;
}

bool p_intersect(const line_t &line, const plane_t &plane, std::pair<std::optional<float>, std::optional<float>> &ptimes) {
    /* if (line.n.dot(plane.normal) > 0 && plane.sided) {
        return false;
    } */
    const float alpha = plane.normal.x * line.n.x + plane.normal.y * line.n.y + plane.normal.z * line.n.z;
    if (alpha == 0) {
        return false;
    }
    ptimes.first = (plane.normal.x * (plane.pos.x - line.pos.x) + 
        plane.normal.y * (plane.pos.y - line.pos.y) +
        plane.normal.z * (plane.pos.z - line.pos.z)) / alpha;
    return true;
}

template <std::uint32_t C>
bool pg_intersect(const line_t &line, const polygon_t<C> &polygon, std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc) {
    static thread_local std::pair<std::optional<float>, std::optional<float>> thisout;
    thisout.first.reset();
    thisout.second.reset();
    const plane_t plane = polygon.plane();
    if (!p_intersect(line, plane, thisout)) {
        return false;
    }
    if (!thisout.first.has_value()) {
        notcurses_stop(nc);
        ERR_EXIT("p_intersect returned true but did not return an intersection point");
    }
    float t = thisout.first.value();
    const vec3_t pos = line.f(t);
    bool b = true;
    for (std::uint32_t i = 0; i < C; i++) {
        const vec3_t a = (polygon.v[i] + polygon.v[(i + 1) % C]) / 2.0f;
        b = b && std::signbit((a - plane.pos).dot(pos - a));
    }
    if (b) {
        ptimes.first = t;
        return true;
    }
    return false;
}

bool cl_intersect(const line_t &line, const cylinder_t &cylinder, std::pair<std::optional<float>, std::optional<float>> &ptimes) {
    const float malpha = cylinder.l.n.x * line.n.x + cylinder.l.n.y * line.n.y;
    const float alpha = malpha * (malpha + 2 * cylinder.l.n.z * line.n.z);
    const float beta = 2 * cylinder.l.n.x * cylinder.l.n.x * line.n.x * (line.pos.x - cylinder.l.pos.x) +
        2 * cylinder.l.n.y * cylinder.l.n.y * line.n.y * (line.pos.y - cylinder.l.pos.y) +
        2 * cylinder.l.n.z * cylinder.l.n.z * line.n.z * (line.pos.z - cylinder.l.pos.z) +
        cylinder.l.n.x * cylinder.l.n.y * (line.pos.x * line.n.y + line.pos.y * line.n.x - cylinder.l.pos.x * line.n.y - cylinder.l.pos.y * line.n.x) + 
        cylinder.l.n.z * cylinder.l.n.y * (line.pos.z * line.n.y + line.pos.y * line.n.z - cylinder.l.pos.z * line.n.y - cylinder.l.pos.y * line.n.z) + 
        cylinder.l.n.z * cylinder.l.n.x * (line.pos.z * line.n.x + line.pos.x * line.n.z - cylinder.l.pos.z * line.n.x - cylinder.l.pos.x * line.n.z);
    const float gamma = cylinder.l.n.x * cylinder.l.n.x * (cylinder.l.pos.x - line.pos.x) * (cylinder.l.pos.x - line.pos.x) +
        cylinder.l.n.y * cylinder.l.n.y * (cylinder.l.pos.y - line.pos.y) * (cylinder.l.pos.y - line.pos.y) +
        cylinder.l.n.z * cylinder.l.n.z * (cylinder.l.pos.z - line.pos.z) * (cylinder.l.pos.z - line.pos.z) +
        cylinder.l.n.x * cylinder.l.n.y;
    return false;
}

bool g_intersect(const line_t &line, const gobj_t &gobj, std::pair<std::optional<float>, std::optional<float>> &ptimes, struct notcurses *nc) {
    bool i = false;
    if (gobj.obj.index() == gtype_sphere) {
        i = i || s_intersect(line, std::get<sphere_t>(gobj.obj), ptimes);
    } else if (gobj.obj.index() == gtype_plane) {
        i = i || p_intersect(line, std::get<plane_t>(gobj.obj), ptimes);
    } else if (gobj.obj.index() == gtype_rect) {
        i = i || pg_intersect(line, std::get<rect_t>(gobj.obj), ptimes, nc);
    } else if (gobj.obj.index() == gtype_triangle) {
        i = i || pg_intersect(line, std::get<triangle_t>(gobj.obj), ptimes, nc);
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
    return i;
}

vec3_t p_reflect(const vec3_t &vec, const plane_t &plane) {
    /* vec3_t effective_normal = (plane.normal * (2 * (plane.normal.dot(vec) < 0) - 1)).normalized(); */
    vec3_t effective_normal = plane.normal;
    return vec - effective_normal * vec.dot(effective_normal) * 2.0f;
}

vec3_t s_reflect(const vec3_t &vec, const sphere_t &sphere, const vec3_t &pos) {
    return p_reflect(vec, sphere.normal_plane(pos));
}

vec3_t g_reflect(const vec3_t &vec, const gobj_t &gobj, const vec3_t &pos, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        return s_reflect(vec, std::get<sphere_t>(gobj.obj), pos);
    } else if (gobj.obj.index() == gtype_plane) {
        return p_reflect(vec, std::get<plane_t>(gobj.obj));
    } else if (gobj.obj.index() == gtype_rect) {
        return p_reflect(vec, std::get<rect_t>(gobj.obj).plane());
    } else if (gobj.obj.index() == gtype_triangle) {
        return p_reflect(vec, std::get<triangle_t>(gobj.obj).plane());
    }
    notcurses_stop(nc);
    ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
}

vec3_t rotated_x_about(const vec3_t &vec, const vec3_t &pos, float theta) {
    return (vec - pos).x_rotated(theta) + pos;
}

vec3_t rotated_y_about(const vec3_t &vec, const vec3_t &pos, float theta) {
    return (vec - pos).y_rotated(theta) + pos;
}


wchar_t scene_t::get_gradient(float x) {
    return gradient.at(std::clamp<std::size_t>(static_cast<std::int32_t>(std::round(x)), 0, gradient.length() - 1));
}

void scene_t::intersect_ray(const line_t &line, const gobj_t *last_obj, gobj_t *&closest_obj, float &closest_t, struct notcurses *nc) {
    static thread_local std::pair<std::optional<float>, std::optional<float>> ptimes;
    /* first, determine the closest intersection point out of all objects on the ray */
    closest_t = std::numeric_limits<float>::max();
    vec3_t closest_pos{closest_t, closest_t, closest_t}; /* closest_t is max value right now */
    /* prgobj = pointer to (potential) reflecting gobj */
    for (gobj_t *prgobj : objects) {
        if (prgobj == last_obj || prgobj->hidden) { continue; }
        ptimes.first.reset();
        ptimes.second.reset();
        if (!g_intersect(line, *prgobj, ptimes, nc)) {
            continue;
        }
        float current_t = optional_min(ptimes, nc);
        /* make sure that it is "forward" on the ray, since light has directprismion */
        if (current_t < 0.0f) { continue; }
        if (current_t < closest_t) {
            closest_t = current_t;
            closest_obj = prgobj;
        }
    }
}

std::uint64_t scene_t::render_ray(const line_t &ray, rgb_t &outcolor, wchar_t &outchar, float init_distance, struct notcurses *nc) {
    std::uint64_t total_light_bounces = 0;
    rgb_t total_color;
    float total_applied_light = 0.0f;
    std::int32_t samples = 0;

    for (samples = 0; samples < samples_per_ray; samples++) {
        std::uint64_t light_bounces = 0;
        const gobj_t *original_obj = nullptr;
        const gobj_t *light = nullptr;
        const gobj_t *last_obj = nullptr;
        rgb_t color{255, 255, 255};
        rgb_t light_color;
        scene_t *current_scene = this;
        line_t line = ray;
        float total_distance = init_distance;

        for (light_bounces = 0; light_bounces < max_light_bounces; light_bounces++) {
            gobj_t *closest_obj = nullptr;
            float closest_t = std::numeric_limits<float>::max();
            current_scene->intersect_ray(line, last_obj, closest_obj, closest_t, nc);
            vec3_t closest_pos = line.f(closest_t);

            if (closest_obj == nullptr) {
                break;
            }
            last_obj = closest_obj;

            if (original_obj == nullptr && !closest_obj->mirror && closest_obj->portal == nullptr) {
                original_obj = closest_obj;
                color = closest_obj->texture(*closest_obj, closest_pos);
            }
            total_distance += (line.pos - closest_pos).mod();
            line.pos = closest_pos;
            if (closest_obj->portal != nullptr) {
                current_scene = closest_obj->portal;
                continue;
            }
            if (closest_obj->transparent) {
                float opacity = closest_obj->opacity_texture(*closest_obj, closest_pos);
                rgb_t tcolor = multiplier(color, opacity) + multiplier(closest_obj->texture(*closest_obj, closest_pos), 1.0f - opacity);
                color.x = static_cast<std::int32_t>(static_cast<float>(color.x * tcolor.x) / 255.0f);
                color.y = static_cast<std::int32_t>(static_cast<float>(color.y * tcolor.y) / 255.0f);
                color.z = static_cast<std::int32_t>(static_cast<float>(color.z * tcolor.z) / 255.0f);
                continue;
            }

            if (closest_obj->light) {
                /* done, we are not reflecting off a light */
                light = closest_obj;
                light_color = light->texture(*light, closest_pos);
                break;
            }
            vec3_t newvec = g_reflect(line.n, *closest_obj, closest_pos, nc);
            float roughness = closest_obj->roughness_texture(*closest_obj, closest_pos);
            newvec = newvec.x_rotated(get_random_real(-roughness / 2.0f, roughness / 2.0f));
            newvec = newvec.y_rotated(get_random_real(-roughness / 2.0f, roughness / 2.0f));
            line.n = newvec;
        }
        total_light_bounces += light_bounces;


        if (last_obj != nullptr) {
            if (light == nullptr && light_bounces > 0) { /* never got illuminated */
                total_applied_light += 1.0f / static_cast<float>(gradient.size());
                total_color += multiplier(color, minimum_color_multiplier);
            } else {
                float applied_light = light->strength(total_distance);
                if (light_bounces > 0) {
                    color.x = static_cast<std::int32_t>(static_cast<float>(color.x * light_color.x) / 255.0f);
                    color.y = static_cast<std::int32_t>(static_cast<float>(color.y * light_color.y) / 255.0f);
                    color.z = static_cast<std::int32_t>(static_cast<float>(color.z * light_color.z) / 255.0f);
                }
                total_applied_light += applied_light;
                total_color += color;
            }
        }
    }
    
    total_applied_light /= static_cast<float>(samples);
    outchar = get_gradient(static_cast<float>(gradient.size() - 1) * total_applied_light);
    total_applied_light /= static_cast<float>(samples);
    outcolor = multiplier(total_color, total_applied_light);

    return total_light_bounces;
}

void add_rectprism_light(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror, const decltype(gobj_t::strength) &strength) {
    for (const rect_t &rect : {rectprism.bottom, rectprism.top, rectprism.left, rectprism.right, rectprism.back, rectprism.front}) {
        scene.objects.push_back(new gobj_t{
            rect,
            color,
            mirror,
            true,
            strength
        });
    }
}

void add_rectprism(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror) {
    for (const rect_t &rect : {rectprism.bottom, rectprism.top, rectprism.left, rectprism.right, rectprism.back, rectprism.front}) {
        scene.objects.push_back(new gobj_t{
            rect,
            color,
            mirror
        });
    }
}
