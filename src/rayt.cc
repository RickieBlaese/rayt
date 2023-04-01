#include <iostream>
#include <utility>
#include <chrono>
#include <numbers>
#include <map>
#include <vector>
#include <algorithm>
#include <cinttypes>
#include <cmath>

#include <ncurses.h>


WINDOW *init_ncurses() {
    WINDOW *win = initscr();
    noecho();
    keypad(win, true);
    cbreak();
    return win;
}

void deinit_ncurses(WINDOW *win) {
    wclear(win);
    endwin();
}


struct vec_t {
    float x = 0.0f, y = 0.0f, z = 0.0f;

    vec_t() = default;
    vec_t(float x, float y, float z) : x(x), y(y), z(z) {}
    vec_t(const vec_t& other) = default;
    ~vec_t() = default;
    bool operator==(const vec_t& other) const { return x == other.x && y == other.y && z == other.z; }
    vec_t& operator=(const vec_t& other) = default;

    vec_t operator+(const vec_t& other) const { return {x + other.x, y + other.y, z + other.z}; }
    vec_t& operator+=(const vec_t& other) { return *this = *this + other; }
    vec_t operator-() const { return {-x, -y, -z}; }
    vec_t operator-(const vec_t& other) const { return *this + (-other); }
    vec_t& operator-=(const vec_t& other) { return *this = *this - other; }
    vec_t operator*(const float v) const { return {x * v, y * v, z * v}; }
    vec_t operator/(const float v) const { return {x / v, y / v, z / v}; }
    float mod() const { return std::sqrt(x * x + y * y + z * z); }
    float dot(const vec_t& other) const { return x * other.x + y * other.y + z * other.z; }
    vec_t normalized() const { return *this / mod(); }
};

/*
vec_t operator"" _vec(const char *s) {
    std::string in = s;
    if (in.find('{') == std::string::npos || in.find('}') == std::string::npos || (in.find('{') > in.find('}'))) {
        throw std::runtime_error("cannot convert \"" + in + "\" to vec_t");
    }
    in = in.substr(in.find('{'), in.find('}') - in.find('}'));
    std::vector<std::int32_t> nums;
    for (decltype(in)::iterator it = in.begin(); it != in.end(); it++) {
        nums.push_back(
    return {std::atoi(in
} */



struct line_t {
    vec_t pos;
    float a = 0, b = 0, c = 0;

    vec_t f(float t) const {
        return {pos.x + a * t, pos.y + b * t, pos.z + c * t};
    }
};

struct sphere_t {
    vec_t pos;
    float r = 0;
};

void intersect(const line_t& line, const sphere_t& sphere, std::vector<vec_t>& out) {
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

line_t between(const vec_t& p1, const vec_t& p2) {
    return {p1, p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
}

std::vector<vec_t>::const_iterator closest_vec(const std::vector<vec_t>& vecs, const vec_t& pos) {
    if (vecs.empty()) {
        deinit_ncurses(stdscr);
        fprintf(stderr, "error: " __FILE__ ":%s():%i: empty vecs passed\n", __func__, __LINE__);
        exit(EXIT_FAILURE);
    }

    float mindist = (pos - vecs[0]).mod();
    auto minvec = vecs.begin();
    for (auto it = vecs.begin() + 1; it != vecs.end(); it++) {
        if ((pos - *it).mod() < mindist) {
            mindist = (pos - *it).mod();
            minvec = it;
        }
    }
    return minvec;
}
    

struct vec4_t {
    vec_t a, b;
};

/* acts like point light */
struct light_t {
    sphere_t sphere;
    float strength = 0.0f;
};

struct obj_t {
    sphere_t sphere;
    float smoothness = 1.0f;
};

struct scene_t {
    std::vector<obj_t> objects;
    std::vector<light_t> lights;
};

struct char_ex_info_t {
    char ch;
    bool is_bold = false;
};


std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}


int main() {
    constexpr float fps = 60.0f;

    const std::string gradient = ".._,'`^\"-~:;=l!i><+?|][}{)(\\/1trxnuvczjfXYUJICLQO0Zmwqpdbkhao*#MW&8%B@$";
    const std::int32_t gradient_length = gradient.length();
    const auto get_gradient = [&gradient, &gradient_length](float x) -> char {
        return gradient[std::clamp<std::int32_t>(static_cast<std::int32_t>(std::round(x)), 0, gradient_length - 1)];
    };

    const float wait_per_frame = 1.0f / fps;
    const std::uint64_t wait_us_per_frame = 1'000'000ULL / fps;

    WINDOW *win = init_ncurses();
    wclear(win);

    vec_t camera_pos(0, 0, 0);

    const auto gaussian_d = [](float x, float y) { return (x * x + y * y) * 100 * std::pow(std::numbers::e, -x * x - y * y); };

    nodelay(win, true);
    curs_set(0);
    std::uint64_t begin_time = get_current_time();
    std::uint64_t last_time  = begin_time;
    const float begin_draw_dist = 30.0f;

    /* aspect ratio ? */
    const float x_mul = 1.0f;
    const float y_mul = 1.0f;

    const float z_move_speed = 20.0f;
    const float x_move_speed = 20.0f;
    const float y_move_speed = 20.0f;

    scene_t scene;
    scene.objects.push_back(obj_t{sphere_t{vec_t(0, 0, 30), 10.0f}, 0.1f});
    scene.lights.push_back(light_t{sphere_t{vec_t(0, 15, 30), 1.0f}, 50.0f});
    const vec_t light_orig_pos = scene.lights[0].sphere.pos;

    std::vector<vec_t> ipoints;
    std::map<float, char_ex_info_t> dist_to_chars;
    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        last_time = get_current_time();
        scene.lights[0].sphere.pos.x = light_orig_pos.x + 10 * std::cos((last_time - begin_time) / 1'000'000.0);
        scene.lights[0].sphere.pos.z = light_orig_pos.z + 10 * std::sin((last_time - begin_time) / 1'000'000.0);
        chtype chin = wgetch(win);
        if (chin == '\t') { break; }
        else if (chin == 'q') { camera_pos.y += wait_per_frame * y_move_speed; }
        else if (chin == 'e') { camera_pos.y -= wait_per_frame * y_move_speed; }
        else if (chin == 'w') { camera_pos.z += wait_per_frame * z_move_speed; }
        else if (chin == 's') { camera_pos.z -= wait_per_frame * z_move_speed; }
        else if (chin == 'd') { camera_pos.x += wait_per_frame * x_move_speed; }
        else if (chin == 'a') { camera_pos.x -= wait_per_frame * x_move_speed; }

        for (std::int32_t i = 0; i < LINES; i++) {
            for (std::int32_t j = 0; j < COLS / 2.0f; j++) {
                /* position for i is flipped since ncurses says y down is positive while we want y up is positive */
                vec_t curpos = vec_t((j - COLS / 4.0f) * x_mul, (-i + LINES / 2.0f) * y_mul, begin_draw_dist);
                line_t ray = between(camera_pos, camera_pos + curpos);
                char out = ' ';
                dist_to_chars.clear();

                /* render spheres with appropriate lighting */
                for (const obj_t& object : scene.objects) {
                    ipoints.clear();
                    intersect(ray, object.sphere, ipoints);
                    if (!ipoints.empty()) {
                        vec_t sphere_minvec = *closest_vec(ipoints, camera_pos);
                        /* nonpermanent and bad solution for camera being in front of sphere */
                        if (sphere_minvec.z > camera_pos.z) {
                            out = '.';
                            /* checking to make sure it is the visible intersection
                             * i.e. the closest intersect with sphere is equal or close
                             * enough to the closest intersect with light source */
                            float applied_light = 0.0f;
                            for (const light_t& light : scene.lights) {
                                ipoints.clear();
                                const float dotp = (object.sphere.pos - sphere_minvec).normalized().dot((object.sphere.pos - light.sphere.pos).normalized());
                                intersect(between(light.sphere.pos, sphere_minvec), object.sphere, ipoints);
                                /* guaranteed to intersect because by above line it goes through both light and sphere */
                                vec_t lminvec;
                                if (ipoints.empty()) {
                                    lminvec = sphere_minvec;
                                } else {
                                    lminvec = *closest_vec(ipoints, light.sphere.pos);
                                }
                                if (dotp > object.smoothness && (lminvec - sphere_minvec).mod() < 0.05f) {
                                    applied_light += light.strength * dotp;
                                }
                            }
                            char outch = '.';
                            if (applied_light > 0.0f) {
                                outch = get_gradient(gradient_length - applied_light);
                            }
                            dist_to_chars[(camera_pos - sphere_minvec).mod()] = char_ex_info_t{outch, false};
                        }
                    }
                }

                /* render light sources */
                for (const light_t& light : scene.lights) {
                    /* checking intersect with light src */
                    ipoints.clear();
                    intersect(ray, light.sphere, ipoints);
                    if (!ipoints.empty()) {
                        /* find closest point of intersect, we do not want back of sphere */
                        const vec_t minvec = *closest_vec(ipoints, camera_pos);

                        /* nonpermanent and bad solution */
                        if (minvec.z > camera_pos.z) {
                            dist_to_chars[(camera_pos - light.sphere.pos).mod()] = char_ex_info_t{get_gradient(light.strength), true};
                        }
                    }
                }

                if (dist_to_chars.empty()) {
                    mvwaddch(win, i, j * 2, ' ');
                    waddch(win, ' ');
                } else {
                    const char_ex_info_t outch = dist_to_chars.begin()->second;
                    if (outch.is_bold) {
                        wattron(win, A_BOLD);
                    }
                    mvwaddch(win, i, j * 2, outch.ch);
                    waddch(win, outch.ch);
                    if (outch.is_bold) {
                        wattroff(win, A_BOLD);
                    }
                }
            }
        }

        mvwprintw(win, 0, 0, "x: %f", camera_pos.x);
        mvwprintw(win, 1, 0, "y: %f", camera_pos.y);
        mvwprintw(win, 2, 0, "z: %f", camera_pos.z);

        wrefresh(win);
    }


    deinit_ncurses(win);
    return 0;
}
