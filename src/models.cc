#include "models.h"


rgb_t multiplier(const rgb_t &original_color, double k) {
    rgb_t outcolor;
    outcolor.x = static_cast<std::int32_t>(static_cast<double>(original_color.x) * k);
    outcolor.y = static_cast<std::int32_t>(static_cast<double>(original_color.y) * k);
    outcolor.z = static_cast<std::int32_t>(static_cast<double>(original_color.z) * k);
    return outcolor;
}

rgb_t color_multiplier(const rgb_t &a, const rgb_t &b) {
    rgb_t outcolor;
    outcolor.x = static_cast<std::int32_t>(static_cast<double>(a.x) * static_cast<double>(b.x) / 255.0);
    outcolor.y = static_cast<std::int32_t>(static_cast<double>(a.y) * static_cast<double>(b.y) / 255.0);
    outcolor.z = static_cast<std::int32_t>(static_cast<double>(a.z) * static_cast<double>(b.z) / 255.0);
    return outcolor;
}


rgb_t average_colors(const rgb_t &a, const rgb_t &b) {
    return multiplier(a + b, 0.5);
}

rgb_t clamp(const rgb_t &color) {
    return {std::clamp(color.x, 0, 255), std::clamp(color.y, 0, 255), std::clamp(color.z, 0, 255)};
}


inline void sort_by_dist(std::vector<vec3_t> &vecs, const vec3_t &pos) {
    std::sort(vecs.begin(), vecs.end(), [&pos](const vec3_t &a, const vec3_t &b) { return (pos - a).mod() < (pos - b).mod(); });
}

vec3_t line_t::f(double t) const {
    return {pos.x + n.x * t, pos.y + n.y * t, pos.z + n.z * t};
}

line_t line_between(const vec3_t &p1, const vec3_t &p2) {
    return {p1, p2 - p1};
}

plane_t sphere_t::normal_plane(const vec3_t &loc) const {
    /* multiplying by r accounts for inverse spheres, making the normal opposite if r is negative */
    return {loc, ((loc - pos)).normalized()};
}

rectprism_t rectprism_t::reversed_wo() const {
    return {bottom.reversed_wo(), top.reversed_wo(),
        left.reversed_wo(), right.reversed_wo(),
        back.reversed_wo(), front.reversed_wo()};
}

vec3_t rectprism_t::get_pos() const {
    vec3_t acc;
    for (const polygon_t &rect : {bottom, top, left, right, back, front}) {
        acc += rect.get_pos();
    }
    return acc / 6.0;
}

rectprism_t create_rectprism(const vec3_t &pos, const vec3_t &size) {
    return {
        /* bottom */
        polygon_t{{
            pos, pos + vec3_t(0, 0, size.z),
            pos + vec3_t(size.x, 0, size.z), pos + vec3_t(size.x, 0, 0)
        }}.reversed_wo(),
        /* top */
        polygon_t{{
            pos + vec3_t(0, size.y, 0), pos + vec3_t(0, size.y, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, size.y, 0)
        }},
        /* left */
        polygon_t{{
            pos, pos + vec3_t(0, 0, size.z),
            pos + vec3_t(0, size.y, size.z), pos + vec3_t(0, size.y, 0)
        }},
        /* right */
        polygon_t{{
            pos + vec3_t(size.x, 0, 0), pos + vec3_t(size.x, 0, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, size.y, 0)
        }}.reversed_wo(),
        /* back */
        polygon_t{{
            pos, pos + vec3_t(0, size.y, 0),
            pos + vec3_t(size.x, size.y, 0), pos + vec3_t(size.x, 0, 0)
        }},
        /* front */
        polygon_t{{
            pos + vec3_t(0, 0, size.z), pos + vec3_t(0, size.y, size.z),
            pos + vec3_t(size.x, size.y, size.z), pos + vec3_t(size.x, 0, size.z)
        }}.reversed_wo()
    };
}

plane_t cylinder_t::normal_plane(const vec3_t &loc) const {
    const double t = - (l.pos - loc).dot(l.n) / l.n.dot(l.n);
    return {loc, (loc - l.f(t)).normalized()};
}

double default_light_strength([[maybe_unused]] double x) {
    /* this is just a smooth step-down function from a to b */
    const double a = 0, b = 200;
    x = std::clamp<double>(x, a, b);
    const double alpha = std::pow(std::numbers::e_v<double>, - (b - a) / (x - a));
    const double beta  = std::pow(std::numbers::e_v<double>, - (b - a) / (b - x));
    return 1.0 - alpha / (alpha + beta);
}

rgb_t default_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.color;
}

double default_opacity_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.opacity;
}

double default_roughness_texture([[maybe_unused]] const gobj_t &self, [[maybe_unused]] const vec3_t &pos) {
    return self.roughness;
}

vec3_t gobj_get_pos(const gobj_t &gobj, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        return std::get<sphere_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_plane) {
        return std::get<plane_t>(gobj.obj).pos;
    } else if (gobj.obj.index() == gtype_polygon) {
        return std::get<polygon_t>(gobj.obj).plane().pos;
    } else if (gobj.obj.index() == gtype_cylinder) {
        return std::get<cylinder_t>(gobj.obj).l.pos;
    }
    notcurses_stop(nc);
    ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
}

void gobj_set_pos(gobj_t &gobj, const vec3_t &pos, struct notcurses *nc) {
    if (gobj.obj.index() == gtype_sphere) {
        std::get<sphere_t>(gobj.obj).pos = pos;
    } else if (gobj.obj.index() == gtype_plane) {
        std::get<plane_t>(gobj.obj).pos = pos;
    } else if (gobj.obj.index() == gtype_polygon) {
        auto &polygon = std::get<polygon_t>(gobj.obj);
        const vec3_t d = pos - polygon.plane().pos;
        for (vec3_t &v : polygon.v) {
            v += d;
        }
    } else if (gobj.obj.index() == gtype_cylinder) {
        std::get<cylinder_t>(gobj.obj).l.pos = pos;
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
}

bool s_intersect(const line_t &line, const sphere_t &sphere, std::pair<std::optional<double>, std::optional<double>> &ptimes) {
    const double alpha = (line.n.x * line.n.x) + (line.n.y * line.n.y) + (line.n.z * line.n.z);
    const double beta  = 2 * (
        line.n.x * (line.pos.x - sphere.pos.x) +
        line.n.y * (line.pos.y - sphere.pos.y) +
        line.n.z * (line.pos.z - sphere.pos.z));
    const double gamma = -2 * ((line.pos.x * sphere.pos.x) + (line.pos.y * sphere.pos.y) + (line.pos.z * sphere.pos.z)) +
        line.pos.x * line.pos.x + sphere.pos.x * sphere.pos.x +
        line.pos.y * line.pos.y + sphere.pos.y * sphere.pos.y +
        line.pos.z * line.pos.z + sphere.pos.z * sphere.pos.z -
        (sphere.r * sphere.r);
    const double discr = (beta * beta) - 4 * alpha * gamma;
    if (discr < 0) {
        return false;
    }
    ptimes.first  = (-beta + std::sqrt(discr)) / (2 * alpha);
    ptimes.second = (-beta - std::sqrt(discr)) / (2 * alpha);
    return true;
}

bool p_intersect(const line_t &line, const plane_t &plane, std::pair<std::optional<double>, std::optional<double>> &ptimes) {
    /* if (line.n.dot(plane.normal) > 0 && plane.sided) {
        return false;
    } */
    const double alpha = plane.normal.x * line.n.x + plane.normal.y * line.n.y + plane.normal.z * line.n.z;
    if (alpha == 0) {
        return false;
    }
    ptimes.first = (plane.normal.x * (plane.pos.x - line.pos.x) + 
        plane.normal.y * (plane.pos.y - line.pos.y) +
        plane.normal.z * (plane.pos.z - line.pos.z)) / alpha;
    return true;
}

bool pg_intersect(const line_t &line, const polygon_t &polygon, std::pair<std::optional<double>, std::optional<double>> &ptimes, struct notcurses *nc) {
    static thread_local std::pair<std::optional<double>, std::optional<double>> thisout;
    thisout.first.reset();
    thisout.second.reset();
    const plane_t plane = polygon.plane();
    /* std::pair<std::optional<double>, std::optional<double>> tout;
    if (!s_intersect(line, sphere_t{plane.pos, (polygon.v[0] - plane.pos).mod()}, tout)) {
        return false;
    } */
    if (!p_intersect(line, plane, thisout)) {
        return false;
    }
    if (!thisout.first.has_value()) {
        notcurses_stop(nc);
        ERR_EXIT("p_intersect returned true but did not return an intersection point");
    }
    double t = thisout.first.value();
    const vec3_t pos = line.f(t);
    bool f = true;
    for (std::uint32_t i = 0; i < polygon.v.size(); i++) {
        const vec3_t a = polygon.v[i];
        const vec3_t b = polygon.v[(i + 1) % polygon.v.size()];
        const vec3_t ab = (a + b) / 2.0;
        f = f && std::signbit((b - a).cross(plane.normal).dot(pos - ab));
    }
    if (f) {
        ptimes.first = t;
        return true;
    }
    return false;
}

bool cl_intersect(const line_t &line, const cylinder_t &cylinder, std::pair<std::optional<double>, std::optional<double>> &ptimes) {
    const double bc = cylinder.l.n.z * line.n.y - line.n.z * cylinder.l.n.y;
    const double ca = cylinder.l.n.x * line.n.z - line.n.x * cylinder.l.n.z;
    const double ab = cylinder.l.n.y * line.n.x - line.n.y * cylinder.l.n.x;

    const double zy = cylinder.l.n.y * (cylinder.l.pos.z - line.pos.z) - cylinder.l.n.z * (cylinder.l.pos.y - line.pos.y);
    const double xz = cylinder.l.n.z * (cylinder.l.pos.x - line.pos.x) - cylinder.l.n.x * (cylinder.l.pos.z - line.pos.z);
    const double yx = cylinder.l.n.x * (cylinder.l.pos.y - line.pos.y) - cylinder.l.n.y * (cylinder.l.pos.x - line.pos.x);
    const double alpha = bc * bc + ca * ca + ab * ab;
    const double beta = 2 * (bc * zy + ca * xz + ab * yx);
    const double gamma = zy * zy + xz * xz + yx * yx - (cylinder.l.n.x * cylinder.l.n.x + cylinder.l.n.y * cylinder.l.n.y + cylinder.l.n.z * cylinder.l.n.z) * cylinder.r * cylinder.r;
    const double discr = (beta * beta) - 4 * alpha * gamma;
    if (discr < 0) {
        return false;
    }
    const double ta = (-beta + std::sqrt(discr)) / (2 * alpha);
    const double tb = (-beta - std::sqrt(discr)) / (2 * alpha);
    const vec3_t tap = line.f(ta);
    const vec3_t tbp = line.f(tb);
    const double cta = - (cylinder.l.pos - tbp).dot(cylinder.l.n) / cylinder.l.n.dot(cylinder.l.n);
    const double ctb = - (cylinder.l.pos - tap).dot(cylinder.l.n) / cylinder.l.n.dot(cylinder.l.n);
    if (cylinder.bounds.has_value()) {
        if (tb > ta && (ctb > cylinder.bounds.value().second || cylinder.bounds.value().first > ctb)) {
            return false;
        }
        if (ta > tb && (cta > cylinder.bounds.value().second || cylinder.bounds.value().first > cta)) {
            return false;
        }
    }
    ptimes.first = ta;
    ptimes.second = tb;
    return true;
}

bool g_intersect(const line_t &line, const gobj_t &gobj, std::pair<std::optional<double>, std::optional<double>> &ptimes, struct notcurses *nc) {
    bool i = false;
    if (gobj.obj.index() == gtype_sphere) {
        i = s_intersect(line, std::get<sphere_t>(gobj.obj), ptimes);
    } else if (gobj.obj.index() == gtype_plane) {
        i = p_intersect(line, std::get<plane_t>(gobj.obj), ptimes);
    } else if (gobj.obj.index() == gtype_polygon) {
        i = pg_intersect(line, std::get<polygon_t>(gobj.obj), ptimes, nc);
    } else if (gobj.obj.index() == gtype_cylinder) {
        i = cl_intersect(line, std::get<cylinder_t>(gobj.obj), ptimes);
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
    return i;
}

vec3_t p_reflect(const vec3_t &vec, const plane_t &plane, vec3_t &normal) {
    /* normal = (plane.normal * (2 * (plane.normal.dot(vec) <= 0) - 1)).normalized(); */
    normal = plane.normal;
    return vec - normal * vec.dot(normal) * 2.0;
}

vec3_t g_reflect(const vec3_t &vec, const gobj_t &gobj, const vec3_t &pos, struct notcurses *nc, vec3_t &normal) {
    vec3_t newvec;
    if (gobj.obj.index() == gtype_sphere) {
        newvec = p_reflect(vec, std::get<sphere_t>(gobj.obj).normal_plane(pos), normal);
    } else if (gobj.obj.index() == gtype_plane) {
        newvec = p_reflect(vec, std::get<plane_t>(gobj.obj), normal);
    } else if (gobj.obj.index() == gtype_polygon) {
        newvec = p_reflect(vec, std::get<polygon_t>(gobj.obj).plane(), normal);
    } else if (gobj.obj.index() == gtype_cylinder) {
        newvec = p_reflect(vec, std::get<cylinder_t>(gobj.obj).normal_plane(pos), normal);
    } else {
        notcurses_stop(nc);
        ERR_EXIT("bad variant index: gobj.obj.index() = %zu", gobj.obj.index());
    }
    double roughness = gobj.roughness_texture(gobj, pos);
    vec3_t rvec(0, 0, 1);
    rvec = rvec.y_rotated(get_random_real(0.0, std::numbers::pi_v<double> * 2.0));
    rvec = rvec.x_rotated(get_random_real(0.0, std::numbers::pi_v<double>));
    rvec = rvec * (-1 * static_cast<bool>(normal.dot(rvec) < 0)); /* if it's not in the right hemisphere, just negate */
    return (newvec * (1.0 - roughness) + rvec * roughness).normalized() * vec.mod();
}

vec3_t rotated_x_about(const vec3_t &vec, const vec3_t &pos, double theta) {
    return (vec - pos).x_rotated(theta) + pos;
}

vec3_t rotated_y_about(const vec3_t &vec, const vec3_t &pos, double theta) {
    return (vec - pos).y_rotated(theta) + pos;
}


wchar_t scene_t::get_gradient(double x) {
    return gradient.at(std::clamp<std::size_t>(static_cast<std::int32_t>(std::round(x)), 0, gradient.length() - 1));
}

void scene_t::intersect_ray(const line_t &line, const gobj_t *last_obj, gobj_t *&closest_obj, double &closest_t, struct notcurses *nc) {
    static thread_local std::pair<std::optional<double>, std::optional<double>> ptimes;
    /* first, determine the closest intersection point out of all objects on the ray */
    closest_t = std::numeric_limits<double>::max();
    vec3_t closest_pos{closest_t, closest_t, closest_t}; /* closest_t is max value right now */
    /* pigobj = pointer to (potential) intersecting gobj */
    for (gobj_t *pigobj : objects) {
        if (pigobj == last_obj || pigobj->hidden) { continue; }
        ptimes.first.reset();
        ptimes.second.reset();
        if (!g_intersect(line, *pigobj, ptimes, nc)) {
            continue;
        }
        double current_t = optional_min(ptimes, nc);
        /* make sure that it is "forward" on the ray, since light has direction */
        if (current_t < 0.0) { continue; }
        if (current_t < closest_t) {
            closest_t = current_t;
            closest_obj = pigobj;
        }
    }
}

std::uint64_t scene_t::render_ray(const line_t &ray, rgb_t &outcolor, wchar_t &outchar, double &applied_light, struct notcurses *nc) {
    std::uint64_t total_light_bounces = 0;
    rgb_t total_color;
    double total_applied_light = 0.0;
    std::int32_t samples = 0;
    static thread_local std::vector<intersection_t> intersections;

    for (samples = 0; samples < samples_per_ray; samples++) {
        std::uint64_t light_bounces = 0;
        const gobj_t *light = nullptr;
        const gobj_t *last_obj = nullptr;
        scene_t *current_scene = this;
        double total_distance = 0.0;
        double total_roughness = 0.0;
        line_t line = ray;
        intersections.clear();

        for (light_bounces = 0; light_bounces < max_light_bounces; light_bounces++) {
            gobj_t *closest_obj = nullptr;
            double closest_t = std::numeric_limits<double>::max();
            current_scene->intersect_ray(line, last_obj, closest_obj, closest_t, nc);
            vec3_t closest_pos = line.f(closest_t);

            if (closest_obj == nullptr) {
                break;
            }
            last_obj = closest_obj;


            total_distance += (line.pos - closest_pos).mod();
            line.pos = closest_pos;

            if (closest_obj->portal != nullptr) {
                current_scene = closest_obj->portal;
                continue;
            }

            intersection_t it;
            it.dist = total_distance;
            it.color = closest_obj->texture(*closest_obj, closest_pos);
            it.transparent = closest_obj->transparent;


            it.opacity = closest_obj->opacity_texture(*closest_obj, closest_pos);
            if (!closest_obj->transparent) {
                vec3_t normal;
                const vec3_t newvec = g_reflect(line.n, *closest_obj, closest_pos, nc, normal);
                it.c = newvec.normalized().dot(normal);
                it.p = 1.0 / (2.0 * std::numbers::pi_v<double>);
                it.brdf = (1.0 - it.roughness) / std::numbers::pi_v<double>;
                it.roughness = closest_obj->roughness_texture(*closest_obj, closest_pos);
                total_roughness += it.roughness;

                line.n = newvec;
            }

            if (closest_obj->light) {
                /* done, we are not reflecting off a light */
                light = closest_obj;
                it.light = true;
                intersections.push_back(it);
                light_bounces++;
                break;
            }
            if (!closest_obj->mirror) {
                intersections.push_back(it);
            }
        }
        total_light_bounces += light_bounces;


        if (last_obj != nullptr) {
            if (!intersections.empty()) {
                if (light == nullptr) { /* never got illuminated */
                    total_color += multiplier(intersections[0].color, minimum_color_multiplier);
                    total_applied_light += minimum_color_multiplier;
                } else if (intersections.rbegin()->light) {
                    rgb_t color = multiplier(intersections.rbegin()->color, light->strength(total_distance));
                    for (auto it = intersections.rbegin() + 1; it != intersections.rend(); it++) {
                        if (it->transparent) {
                            color = multiplier(it->color, it->opacity) + multiplier(color, 1.0 - it->opacity);
                        } else {
                            color = it->color + multiplier(color, it->brdf * it->c / it->p);
                        }
                    }
                    total_applied_light += light->strength(total_distance);
                    total_color += color;
                }
            }
        }
        /* if we never encountered anything that would change based on roughness then don't waste more samples */
        if (total_roughness == 0.0 || last_obj == nullptr) {
            samples++;
            break;
        }
    }
    
    if (total_applied_light > 0.0) {
        total_applied_light /= static_cast<double>(samples);
         total_color = multiplier(total_color, 1.0 / static_cast<double>(samples));
        outchar = get_gradient(static_cast<double>(gradient.size() - 1) * total_applied_light);
        outcolor = total_color;
        applied_light = total_applied_light;
    } else {
        outchar = L' ';
    }

    return total_light_bounces;
}

void add_rectprism_light(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror, const decltype(gobj_t::strength) &strength) {
    for (const polygon_t &rect : {rectprism.bottom, rectprism.top, rectprism.left, rectprism.right, rectprism.back, rectprism.front}) {
        scene.objects.push_back(new gobj_t{
            .obj = rect,
            .color = color,
            .mirror = mirror,
            .light = true,
            .strength = strength
        });
    }
}

void add_rectprism(scene_t &scene, const rectprism_t &rectprism, const rgb_t &color, bool mirror, double roughness) {
    for (const polygon_t &rect : {rectprism.bottom, rectprism.top, rectprism.left, rectprism.right, rectprism.back, rectprism.front}) {
        scene.objects.push_back(new gobj_t{
            .obj = rect,
            .color = color,
            .mirror = mirror,
            .roughness = roughness
        });
    }
}
