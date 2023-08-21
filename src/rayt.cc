#include <algorithm>
#include <atomic>
#include <barrier>
#include <chrono>
#include <codecvt>
#include <concepts>
#include <fstream>
#include <iostream>
#include <thread>
#include <limits>
#include <locale>
#include <mutex>
#include <random>
#include <sstream>

#include "nlohmann/json.hpp"
#include "models.h"


std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

bool get_input_blocking(std::vector<std::string> &rhistory, const std::string &starting_str, struct ncplane *plane, struct notcurses *nc) {
    std::uint32_t dimy = 0, dimx = 0;
    ncplane_dim_yx(plane, &dimy, &dimx);

    std::vector<std::string> history = rhistory;
    history.emplace_back();
    std::uint32_t cpos = 0;
    std::uint32_t history_pos = history.size() - 1;
#define crout() history[history_pos]
    ncplane_options opts;
    opts.flags = NCPLANE_OPTION_HORALIGNED;
    opts.x = 1;
    opts.y = static_cast<std::int32_t>(dimy) - 1;
    opts.rows = 1;
    opts.cols = dimx;
    struct ncplane *tplane = ncplane_create(plane, &opts);
    ncplane_putstr_yx(tplane, 0, 0, starting_str.c_str());
    ncplane_putchar(tplane, ' ');
    notcurses_cursor_enable(nc, static_cast<std::int32_t>(dimy) - 1, 1);
    notcurses_render(nc);
    bool broke = false;
    while (true) {
        /* std::string thiscrout() = starting_str + crout() + std::string(" ", max_length - ccount); */
        ncplane_erase_region(tplane, 0, static_cast<std::int32_t>(starting_str.length()), 1, static_cast<std::int32_t>(dimx - starting_str.length()));
        ncplane_putstr_yx(tplane, 0, static_cast<std::int32_t>(starting_str.length()), crout().c_str());
        notcurses_cursor_enable(nc, static_cast<std::int32_t>(dimy) - 1, static_cast<std::int32_t>(cpos + starting_str.length()));
        ncplane_putchar_yx(tplane, 0, static_cast<std::int32_t>(starting_str.length() + crout().length()), ' ');
        notcurses_render(nc);
        ncinput ni;
        std::uint32_t chin = notcurses_get_blocking(nc, &ni);
        if (ni.alt || ni.ctrl) {
            continue;
        }
        if (chin == NCKEY_BACKSPACE) {
            if (cpos > 0) {
                cpos--;
                crout().erase(crout().begin() + cpos);
            }
            continue;
        } else if (chin == NCKEY_DEL) {
            if (cpos < crout().length()) {
                crout().erase(crout().begin() + cpos);
            }
            continue;
        } else if (chin == NCKEY_ENTER) {
            break;
        } else if (chin == NCKEY_ESC) {
            broke = true;
            break;
        } else if (chin == NCKEY_LEFT) {
            if (cpos > 0) {
                cpos--;
            }
            continue;
        } else if (chin == NCKEY_RIGHT) {
            if (cpos < crout().length()) {
                cpos++;
            }
            continue;
        } else if (chin == NCKEY_UP) {
            if (history_pos > 0) {
                history_pos--;
            }
            cpos = crout().length();
            continue;
        } else if (chin == NCKEY_DOWN) {
            if (history_pos < history.size() - 1) {
                history_pos++;
            }
            cpos = crout().length();
            continue;
        } else if (!isprint(static_cast<char>(chin))) {
            continue;
        }

        crout().insert(crout().begin() + cpos, static_cast<char>(chin));
        cpos++;
    }
    ncplane_destroy(tplane);
    notcurses_cursor_disable(nc);
    rhistory.push_back(crout());
    return broke;
}

void from_json(const nlohmann::json &j, vec3_t &vec) {
    j.at(0).get_to(vec.x);
    j.at(1).get_to(vec.y);
    j.at(2).get_to(vec.z);
}

void from_json(const nlohmann::json &j, rgb_t &color) {
    j.at(0).get_to(color.x);
    j.at(1).get_to(color.y);
    j.at(2).get_to(color.z);
}

void from_json(const nlohmann::json &j, gobj_t &gobj) {
    gobj = gobj_t{
        .color = rgb_t(255, 255, 255),
        .mirror = false,
        .light = false,
        .roughness = 1.0,
        .transparent = false,
        .opacity = 1.0,
        .hidden = false,
    };
    if (j.contains("color")) {
        j.at("color").get_to(gobj.color);
    }
    if (j.contains("mirror")) {
        j.at("mirror").get_to(gobj.mirror);
    }
    if (j.contains("light")) {
        j.at("light").get_to(gobj.light);
    }
    if (j.contains("roughness")) {
        j.at("roughness").get_to(gobj.roughness);
    }
    if (j.contains("transparent")) {
        j.at("transparent").get_to(gobj.transparent);
    }
    if (j.contains("opacity")) {
        j.at("opacity").get_to(gobj.opacity);
    }
    if (j.contains("hidden")) {
        j.at("hidden").get_to(gobj.hidden);
    }
}


int main(int argc, char **argv) {
    /* default values */
    double fps = 120.0;
    std::uint64_t samples_per_ray = 1;
    std::uint32_t max_light_bounces = 20;
    std::int64_t thread_count = 2;

    for (std::int32_t i = 1; i < argc; i++) {
        if (strlen(argv[i]) > 2) {
            if (argv[i][0] != '-') { continue; }
            if (argv[i][1] == 'j') {
                thread_count = std::strtol(&argv[i][2], nullptr, 0);
            } else if (argv[i][1] == 'f') {
                fps = std::strtof(&argv[i][2], nullptr);
            } else if (argv[i][1] == 's') {
                samples_per_ray = std::strtol(&argv[i][2], nullptr, 0);
            } else if (argv[i][1] == 'm') {
                max_light_bounces = std::strtol(&argv[i][2], nullptr, 0);
            }
        }
    }

    std::wstring gradient;
    std::ifstream gradient_file("gradient.txt");
    if (!gradient_file.is_open()) {
        ERR_EXIT("couldn't open gradient.txt\n");
    }
    std::stringstream ss;
    ss << gradient_file.rdbuf();

    /* this part is hacky */
    std::string mbgradient;
    mbgradient = ss.str();
    mbgradient.erase(std::strcspn(mbgradient.c_str(), "\r\n"), mbgradient.size()); /* only count up to the first newline if found, otherwise up to \0 */
    gradient = std::wstring_convert<std::codecvt_utf8<wchar_t>>().from_bytes(mbgradient);

    gradient_file.close();
    if (gradient.empty()) {
        ERR_EXIT("gradient file was empty");
    }

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
    /* const auto gaussian = [](double x, double y) { return std::pow(std::numbers::e, -x * x - y * y); }; */

    std::uint64_t begin_time = get_current_time();
    std::uint64_t last_time  = begin_time;
    double begin_draw_dist = 50.0;

    /* aspect ratio ? */
    const double x_mul = 1.0;
    const double y_mul = 1.0;

    const double z_move_speed = 1.0;
    const double x_move_speed = 1.0;
    const double y_move_speed = 1.0;
    const double shift_multiplier = 5.0;
    const double angle_increment = std::numbers::pi / 50.0;



    world_t world;

    /* --- init scene --- */
    auto *pscene = new scene_t;
    scene_t &scene = *pscene;
    scene.samples_per_ray = samples_per_ray;
    scene.max_light_bounces = max_light_bounces;
    world.scenes.push_back(pscene);
    scene.gradient = gradient;
    scene.camera_ray = line_t{vec3_t(0, 0, 0), vec3_t(0, 0, begin_draw_dist)};
    vec3_t &camera_pos = scene.camera_ray.pos;
    vec3_t &camera_n   = scene.camera_ray.n;
    vec3_t base_camera_n = camera_n;
    gobj_t *grabbed_obj = nullptr;



    const auto light_strength_func = []([[maybe_unused]] double x) {
        return std::pow(std::numbers::e_v<double>, -x/5000.0);
    };


    /* rectprism_t p = create_rectprism(vec3_t(-5, -5, -5), vec3_t(10, 10, 10)).reversed_wo();
    scene.objects.push_back(new gobj_t{
        .obj = p.top,
        .color = {255, 255, 255},
        .light = true,
        .strength = light_strength_func
    });
    scene.objects.push_back(new gobj_t{
        .obj = p.front,
        .color = {255, 255, 255},
        .roughness = 1.0,
    });
    scene.objects.push_back(new gobj_t{
        .obj = p.bottom,
        .color = {255, 255, 255},
        .roughness = 1.0,
    });
    scene.objects.push_back(new gobj_t{
        .obj = p.right,
        .color = {0, 255, 0},
        .roughness = 1.0,
    });
    scene.objects.push_back(new gobj_t{
        .obj = p.left,
        .color = {255, 0, 0},
        .roughness = 1.0,
    });
    scene.objects.push_back(new gobj_t{
        .obj = sphere_t{vec3_t(0, 0, 0), 3},
        .mirror = true
    }); */

    /* scene.objects.push_back(new gobj_t{
        .obj = polygon_t{{vec3_t(-1.0, 4.9, -1.0), vec3_t(1.0, 4.9, -1.0), vec3_t(1.0, 4.9, 1.0), vec3_t(-1.0, 4.9, 1.0)}},
        .color = {255, 255, 255},
        .light = true,
        .strength = light_strength_func
    }); */
    /* scene.objects.push_back(new gobj_t{
        .obj = cylinder_t{line_between(vec3_t(0, 0, 5), vec3_t(1, 0, 5)), 1.8},
        .color = {255, 0, 0},
        .roughness = 1.0
    }); */


    /* scene.objects.push_back(new gobj_t{
        .obj = triangle_t{{vec3_t(0, 0, 0), vec3_t(0, 1, 0), vec3_t(0, 1, 1)}},
        .color = {255, 0, 0},
        .light = true
    }); */

    /* scene.objects.push_back(new gobj_t{
        .obj = sphere_t{vec3_t(-2, -0.5, 10), 2.5},
        .color = {255, 0, 0},
        .roughness = 1.0,
    });
    scene.objects.push_back(new gobj_t{
        .obj = sphere_t{vec3_t(4, -0.5, 10), 2.5},
        .color = {0, 0, 255},
        .roughness = 1.0,
    }); */
    /* scene.objects.push_back(new gobj_t{
        .obj = sphere_t{vec3_t(-5, -0.5, 5), 2.5},
        .color = {0, 0, 255},
        .transparent = true,
        .opacity = 0.5
    }); */
    /* auto portal_scene = new scene_t{
        .objects = {
            new gobj_t{
                .obj = sphere_t{vec3_t(15, 0, 0), 3.0},
                .color = {0, 255, 0}
            },
            new gobj_t{
                .obj = plane_t{vec3_t(0, -3, 0), vec3_t(0, 1, 0)},
                .color = {255, 255, 255}
            },
            new gobj_t{
                .obj = plane_t{vec3_t(0, 50, 0), vec3_t(0, -1, 0)},
                .color = {255, 255, 255},
                .light = true,
                .strength = default_light_strength
            }
        },
        .gradient = gradient
    };
    world.scenes.push_back(portal_scene);
    scene.objects.push_back(new gobj_t{
        .obj = rect_t{{vec3_t(10, -2, -5), vec3_t(10, -2, 5), vec3_t(10, 8, 5),
    vec3_t(10, 8, -5)}}, .portal = portal_scene
    }); */

    scene.objects.push_back(new gobj_t{
        .obj = plane_t{vec3_t(-50, 0, 0), vec3_t(1, 0, 0)},
        .color = {255, 255, 255},
        .light = true,
        .strength = light_strength_func
    });
    scene.objects.push_back(new gobj_t{
        .obj = plane_t{vec3_t(50, 0, 0), vec3_t(-1, 0, 0)},
        .color = {255, 255, 255},
        .light = true,
        .strength = light_strength_func
    });


    std::mutex ncplane_mutex;
    std::barrier render_complete(thread_count + 1);
    std::vector<std::pair<std::int32_t, std::int32_t>> regions; /* y coords */
    std::vector<std::thread*> vthreads;
    vthreads.reserve(thread_count);
    std::vector<double> itimes, thisout;
    std::uint32_t dimy = 0, dimx = 0;
    double camera_rotation_x = 0.0; /* accounts for both x and z */
    double camera_rotation_y = 0.0;
    std::atomic_uint64_t total_bounce_count = 0;
    double grab_ray_time = 0.5;
    std::atomic_bool stop = false;

    std::vector<vec3_t> spoints;

    
    for (std::int32_t i = 0; i < thread_count; i++) {
        auto render_region = [&](const std::atomic_bool *stopt, std::int32_t region_ind) {
            while (!stopt->load()) {
                render_complete.arrive_and_wait();
                const auto [ay, by] = regions[region_ind];
                for (std::int32_t i = ay; i < by; i++) {
                    for (std::int32_t j = 0; j < dimx / 2; j++) {
                        /* position for i is flipped since notcurses says y down is positive while we want y up is positive */
                        vec3_t curpos = vec3_t((j - dimx / 4.0) * x_mul, (-i + dimy / 2.0) * y_mul, 0); /* NOLINT */
                        curpos = curpos.y_rotated(camera_rotation_y).x_rotated(camera_rotation_x); /* is still "at the origin" */
                        line_t ray = line_between(camera_pos, camera_n + curpos); /* t positive is "forward" */
                        wchar_t current_char = L' ';
                        rgb_t current_color;
                        double applied_light = 0.0;

                        total_bounce_count += scene.render_ray(ray, current_color, current_char, applied_light, nc);
                        current_color = clamp(current_color);

                        while (!ncplane_mutex.try_lock()) {;}
                        ncplane_set_fg_rgb8(plane, current_color.x, current_color.y, current_color.z);
                        ncplane_putwc_yx(plane, i, j * 2, current_char);
                        ncplane_putwc(plane, current_char);
                        ncplane_mutex.unlock();
                    }
                }
                render_complete.arrive_and_wait();
            }
        };
        vthreads.push_back(new std::thread(render_region, &stop, i));
    }

    std::random_device device{};
    std::default_random_engine engine(device());

    static std::uniform_real_distribution<double> thetadist(0, 2 * std::numbers::pi);
    static std::uniform_real_distribution<double> phidist(-std::numbers::pi / 2.0, std::numbers::pi / 2.0);
    static std::uniform_int_distribution<std::int32_t> colordist(0, 255);
    std::string last_msg = "";
    bool msg_is_error = false;

    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        std::uint64_t render_time = get_current_time() - last_time;
        last_time = get_current_time();

        while (!ncplane_mutex.try_lock()) {;}
        ncplane_dim_yx(plane, &dimy, &dimx);
        std::uint32_t chin = notcurses_get_nblock(nc, nullptr);
        ncplane_mutex.unlock();

        vec3_t move_d;
        if (chin == NCKEY_ESC) { break; }
        else if (chin == L'q') { move_d.y =  y_move_speed; }
        else if (chin == L'Q') { move_d.y =  y_move_speed * shift_multiplier; }
        else if (chin == L'e') { move_d.y = -y_move_speed; }
        else if (chin == L'E') { move_d.y = -y_move_speed * shift_multiplier; }
        else if (chin == L'w') { move_d.z =  z_move_speed; }
        else if (chin == L'W') { move_d.z =  z_move_speed * shift_multiplier; }
        else if (chin == L's') { move_d.z = -z_move_speed; }
        else if (chin == L'S') { move_d.z = -z_move_speed * shift_multiplier; }
        else if (chin == L'd') { move_d.x =  x_move_speed; }
        else if (chin == L'D') { move_d.x =  x_move_speed * shift_multiplier; }
        else if (chin == L'a') { move_d.x = -x_move_speed; }
        else if (chin == L'A') { move_d.x = -x_move_speed * shift_multiplier; }

        else if (chin == L'h') {
            // random graph generator
            vec3_t newpoint;
            for (const vec3_t &p : spoints) {
                newpoint += p;
            }
            if (spoints.size() > 0) {
                newpoint = newpoint / spoints.size();
            }
            double maxmod = 0.0;
            for (const vec3_t &p : spoints) {
                if (p.mod() > maxmod) {
                    maxmod = p.mod();
                }
            }
            newpoint += vec3_t(1.0, 0.0, 0.0).x_rotated(thetadist(engine)).y_rotated(phidist(engine)) * (maxmod + std::uniform_real_distribution<double>(1.0, 4.0)(engine));
            scene.objects.push_back(new gobj_t{
                .obj = sphere_t{
                    .pos = newpoint,
                    .r = 1
                },
                .color = {colordist(engine), colordist(engine), colordist(engine)},
                .roughness = 1.0
            });
            std::shuffle(spoints.begin(), spoints.end(), engine);
            std::uint32_t joincount = std::uniform_int_distribution<std::uint32_t>(0, spoints.size())(engine);
            for (std::uint32_t i = 0; i < joincount; i++) {
                scene.objects.push_back(new gobj_t{
                    .obj = cylinder_t{
                        .l = line_t{
                            .pos = newpoint,
                            .n = newpoint - spoints[i]
                        },
                        .r = 0.6,
                        .bounds = std::pair<double, double>(-1.0, 0.0)
                    },
                    .color = {colordist(engine), colordist(engine), colordist(engine)},
                    .roughness = 1.0
                });
            }
            spoints.push_back(newpoint);


            // simplex generator
            /* vec3_t newpoint;
            for (const vec3_t &p : spoints) {
                newpoint += p;
            }
            if (spoints.size() > 0) {
                newpoint = newpoint / spoints.size();
            }
            double maxmod = 0.0;
            for (const vec3_t &p : spoints) {
                if (p.mod() > maxmod) {
                    maxmod = p.mod();
                }
            }
            newpoint += vec3_t(1.0, 0.0, 0.0).x_rotated(thetadist(engine)).y_rotated(phidist(engine)) * (maxmod + 1.0);
            spoints.push_back(newpoint);
            for (const vec3_t &p : spoints) {
                scene.objects.push_back(new gobj_t{
                    .obj = cylinder_t{
                        .l = line_t{
                            .pos = newpoint,
                            .n = newpoint - p
                        },
                        .r = 0.2,
                        .bounds = std::pair<double, double>(-1.0, 0.0)
                    },
                    .color = {colordist(engine), colordist(engine), colordist(engine)},
                    .roughness = 1.0
                });
            } */
        }

        else if (chin == L'j' || chin == L'J') { camera_rotation_x += angle_increment; }
        else if (chin == L'l' || chin == L'L') { camera_rotation_x -= angle_increment; }
        else if ((chin == L'i' || chin == L'I') && camera_rotation_y < std::numbers::pi / 2.0 - angle_increment) {
            camera_rotation_y += angle_increment;
        } else if ((chin == L'k' || chin == L'K') && camera_rotation_y > - (std::numbers::pi / 2.0 - angle_increment)) {
            camera_rotation_y -= angle_increment;
        }

        else if (chin == L' ') {
            gobj_t *closest_obj = nullptr;
            double closest_t = std::numeric_limits<double>::max();
            itimes.clear();
            thisout.clear();
            scene.intersect_ray(line_between(camera_pos, camera_n), nullptr, closest_obj, closest_t, nc);
            if (closest_obj != nullptr) {
                delete closest_obj; /* does not modify closest_obj */
                std::erase_if(scene.objects, [&closest_obj](const gobj_t * const &pgobj) { return pgobj == closest_obj; });
            }
        }

        else if (chin == L'c' || chin == L'C') {
            if (grabbed_obj == nullptr) {
                double closest_t = std::numeric_limits<double>::max();
                itimes.clear();
                thisout.clear();
                scene.intersect_ray(line_between(camera_pos, camera_n), nullptr, grabbed_obj, closest_t, nc);
            } else {
                grabbed_obj = nullptr;
            }
        }

        else if (chin == L'x' || chin == L'X') {
            grab_ray_time -= 0.02;
        }
        else if (chin == L'v' || chin == L'V') {
            grab_ray_time += 0.02;
        }

        // all commands
        else if (chin == L';') {
            static thread_local std::vector<std::string> cmdhistory;
            if (get_input_blocking(cmdhistory, ";", plane, nc)) {
                goto end_input;
            }
            std::string in = *(cmdhistory.end() - 1);

            /* ltrim copied from https://stackoverflow.com/questions/216823/how-to-trim-an-stdstring */
            in.erase(in.begin(), std::find_if(in.begin(), in.end(), [](unsigned char ch) { return !std::isspace(ch); }));

            std::size_t spaceloc = in.find_first_of(' ');
            std::string cmd, args;
            if (spaceloc == std::string::npos) {
                spaceloc = in.length(); /* bad but it's what we want */
            } else {
                args = in.substr(spaceloc);
            }
            cmd = in.substr(0, spaceloc);

            /* ltrim */
            args.erase(args.begin(), std::find_if(args.begin(), args.end(), [](unsigned char ch) { return !std::isspace(ch); }));

            if (cmd.empty()) {
                goto end_input;
            }

            if (cmd == "add") {
                nlohmann::json j;
                gobj_t gobj;
                try {
                    j = nlohmann::json::parse(args);
                    gobj = j[0].get<gobj_t>();
                } catch (std::exception &e) {
                    last_msg = std::string("failed to parse json: ") + e.what();
                    msg_is_error = true;
                    goto end_input;
                }
                last_msg = "created gobj";
                msg_is_error = false;
            } else if (cmd == "q") {
                break;
            } else {
                last_msg = "command not recognized";
                    msg_is_error = true;
                goto end_input;
            }
        }


        end_input:


        camera_n = rotated_x_about(rotated_y_about(base_camera_n, camera_pos, camera_rotation_y), camera_pos, camera_rotation_x);

        move_d = move_d.x_rotated(camera_rotation_x);
        camera_pos += move_d;
        camera_n += move_d;
        base_camera_n += move_d;
        std::uint64_t world_objects_count = 0;
        for (const scene_t * const scene : world.scenes) {
            world_objects_count += scene->objects.size();
        }

        if (grabbed_obj != nullptr) {
            gobj_set_pos(*grabbed_obj, line_between(camera_pos, camera_n).f(grab_ray_time), nc);
        }

        regions.clear();
        partition(0, static_cast<std::int32_t>(dimy), static_cast<std::int32_t>(thread_count), regions);

        /* begin rendering */
        total_bounce_count = 0;
        render_complete.arrive_and_wait();

        /* end rendering */
        render_complete.arrive_and_wait();

        while (!ncplane_mutex.try_lock()) {;}
        ncplane_set_fg_rgb8(plane, 255, 255, 255);
        ncplane_printf_yx(plane, 0, 0, "x:  %-+4.2f ", camera_pos.x);
        ncplane_printf_yx(plane, 1, 0, "y:  %-+4.2f ", camera_pos.y);
        ncplane_printf_yx(plane, 2, 0, "z:  %-+4.2f ", camera_pos.z);
        ncplane_printf_yx(plane, 3, 0, "nx: %-+4.2f ", camera_n.x - camera_pos.x);
        ncplane_printf_yx(plane, 4, 0, "ny: %-+4.2f ", camera_n.y - camera_pos.y);
        ncplane_printf_yx(plane, 5, 0, "nz: %-+4.2f ", camera_n.z - camera_pos.z);
        ncplane_printf_yx(plane, 6, 0, "scene obj: %zu ", scene.objects.size());
        ncplane_printf_yx(plane, 7, 0, "world obj: %zu ", world_objects_count);
        ncplane_printf_yx(plane, 8, 0, "render: %lu µs ", render_time);
        ncplane_printf_yx(plane, 9, 0, "fps: %5.4Lf ", 1'000'000.0L / render_time);
        ncplane_printf_yx(plane, 10, 0, "bounces: %lu ", total_bounce_count.load());
        ncplane_printf_yx(plane, 11, 0, "grabbed obj: %p ", grabbed_obj);
        ncplane_putstr_yx(plane, 12, 0, "y regions: ");
        for (std::int32_t i = 0; i < regions.size(); i++) {
            ncplane_printf_yx(plane, i + 13, 0, "%i : %i ", regions[i].first, regions[i].second);
        }
        if (!last_msg.empty()) {
            ncplane_putstr_yx(plane, static_cast<std::int32_t>(dimy) - 2, 0, "> ");
            if (msg_is_error) {
                ncplane_set_fg_rgb8(plane, 255, 10, 70);
            } else {
                ncplane_set_fg_rgb8(plane, 10, 255, 70);
            }
            ncplane_putstr(plane, last_msg.c_str());
            ncplane_set_fg_rgb8(plane, 255, 255, 255);
        }

        notcurses_render(nc);
        ncplane_mutex.unlock();
    }


    /* stopping jthreads */
    stop = true;

    render_complete.arrive_and_wait();
    render_complete.arrive_and_wait();
    for (std::thread *t : vthreads) {
        t->join();
        delete t;
    }


    for (scene_t *scene : world.scenes) {
        for (gobj_t *gobj : scene->objects) {
            delete gobj;
        }
        delete scene;
    }

    bool s = notcurses_cantruecolor(nc);
    ncplane_destroy(plane);
    notcurses_stop(nc);

    std::cout << "cantruecolor: " << std::boolalpha << s << '\n';
    return 0;
}
