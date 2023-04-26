#include <barrier>
#include <chrono>
#include <fstream>
#include <iostream>
#include <notcurses/notcurses.h>
#include <thread>
#include <limits>
#include <mutex>
#include <random>
#include <sstream>

#include <cstdlib>

#include "models.h"



std::uint64_t get_current_time() {
    return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}


int main(int argc, char **argv) {
    float fps = 120.0f;
    std::uint64_t samples_per_ray = 5;
    
    std::int64_t thread_count = 8;
    for (std::int32_t i = 1; i < argc; i++) {
        if (strlen(argv[i]) > 2) {
            if (argv[i][0] != '-') { continue; }
            if (argv[i][1] == 'j') {
                thread_count = std::strtol(&argv[i][2], nullptr, 0);
            } else if (argv[i][1] == 'f') {
                fps = std::strtof(&argv[i][2], nullptr);
            } else if (argv[i][1] == 's') {
                samples_per_ray = std::strtol(&argv[i][2], nullptr, 0);
            }
        }
    }

    /* const std::wstring gradient = L" ._,'`^\"-~:;=!><+?|][}{)(\\/trxnuovczmwaihqpdbkfjl1XYFGHNUJICLQO0Z#MW&8%B@$"; */
    std::wstring gradient;
    std::wifstream gradient_file("gradient.txt");
    std::wstringstream ss;
    ss << gradient_file.rdbuf();
    gradient = ss.str();
    gradient_file.close();

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
    float begin_draw_dist = 50.0f;

    /* aspect ratio ? */
    const float x_mul = 1.0f;
    const float y_mul = 1.0f;

    const float z_move_speed = 1.0f;
    const float x_move_speed = 1.0f;
    const float y_move_speed = 1.0f;
    const float shift_multiplier = 5.0f;
    const float angle_increment = std::numbers::pi / 30.0f;


    world_t world;

    /* --- add objects to scene --- */
    auto *pscene = new scene_t;
    scene_t &scene = *pscene;
    scene.samples_per_ray = samples_per_ray;
    world.scenes.push_back(pscene);
    scene.gradient = gradient;
    scene.camera_ray = line_t{vec3_t(0, 0, 0), vec3_t(0, 0, begin_draw_dist)};
    vec3_t &camera_pos = scene.camera_ray.pos;
    vec3_t &camera_n   = scene.camera_ray.n;
    vec3_t base_camera_n = camera_n;

    
    const auto light_strength_func = []([[maybe_unused]] float x) {
        return std::pow(std::numbers::e_v<float>, -x/3000.0f);
    };

    /* test spheres */
    scene.objects.push_back(new gobj_t{
        sphere_t{vec3_t(0, 0, 30), 20.0f},
        {0, 255, 0}
    });
    scene.objects.push_back(new gobj_t{
        .obj = sphere_t{vec3_t(0, 0, 80), 20.0f},
        .color = {0, 0, 255},
        .transparent = true,
        .opacity = 1.0f
    });


    /* planar mirrors */
    /* scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(0, 0, 105), vec3_t(0, 0, -1)},
            -200, -20, 400, 400
        },
        {0, 0, 0},
        true,
    }); */

    scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(0, -20, 0), vec3_t(0, 1, 0)},
            -200, -200, 400, 400
        },
        {200, 255, 200},
    });
    /* right wall */
    /* scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(90, 0, 0), vec3_t(-1, 0, 0)},
            -400, 0, 800, 800
        },
        {255, 255, 255}
    }); */
    /* back wall */
    /* scene.objects.push_back(new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(0, 0, -90), vec3_t(0, 0, 1)},
            -400, 0, 800, 800
        },
        {255, 255, 255}
    }); */

    /* add_rect(scene, create_rect(true, vec3_t(-20, -20, -50), vec3_t(40, 40, 40)), rgb_t(255, 99, 255), false); */
    auto *portal_scene = new scene_t{
        {
            new gobj_t{sphere_t{vec3_t(-110, 0, 0), 5.0f}, rgb_t(255, 0, 0)},
            new gobj_t{plane_t{vec3_t(0, 100, 0), vec3_t(0, -1, 0)}, rgb_t(255, 255, 255), false, true, light_strength_func},
            new gobj_t{plane_t{vec3_t(0, -5.0f, 0), vec3_t(0, 1, 0)}, rgb_t(255, 255, 255)}
        }
    };
    portal_scene->samples_per_ray = 2.0f;
    world.scenes.push_back(portal_scene);
    portal_scene->gradient = gradient;
    auto *portal_plane = new gobj_t{
        bounded_plane_t{
            plane_t{vec3_t(-100, 0, 0), vec3_t(1, 0, 0)},
            -10, -10, 20, 20
        },
        {255, 255, 255}
    };
    portal_plane->portal = portal_scene;
    scene.objects.push_back(portal_plane);

    /* lights */
    /* add_rect_light(scene, create_rect(false, vec3_t(-800, -800, -800), vec3_t(1600, 1600, 1600)), rgb_t(255, 255, 255), false, light_strength_func); */
    /* scene.objects.push_back(new gobj_t{
        plane_t{vec3_t(0, 500, 0), vec3_t(0, -1, 0)},
        {255, 255, 200},
        false,
        true,
        light_strength_func
    }); */
    scene.objects.push_back(new gobj_t{
        plane_t{vec3_t(0, 400, 0), vec3_t(0, -1, 0)},
        {255, 255, 220},
        false,
        true,
        light_strength_func
    });


    std::mutex ncplane_mutex;
    std::barrier render_complete(thread_count + 1);
    std::vector<std::pair<std::int32_t, std::int32_t>> regions; /* y coords */
    std::vector<std::jthread*> vthreads;
    vthreads.reserve(thread_count);
    std::vector<float> itimes, thisout;
    std::uint32_t dimy = 0, dimx = 0;
    float camera_rotation_x = 0.0f; /* accounts for both x and z */
    float camera_rotation_y = 0.0f;
    std::uint64_t total_bounce_count = 0;

    
    for (std::int32_t i = 0; i < thread_count; i++) {
        auto render_region = [&](const std::stop_token &stoken, std::int32_t region_ind) {
            while (!stoken.stop_requested()) {
                render_complete.arrive_and_wait();
                auto [ay, by] = regions[region_ind];
                for (std::int32_t i = ay; i < by; i++) {
                    for (std::int32_t j = 0; j < dimx / 2; j++) {
                        /* position for i is flipped since notcurses says y down is positive while we want y up is positive */
                        vec3_t curpos = vec3_t((j - dimx / 4.0f) * x_mul, (-i + dimy / 2.0f) * y_mul, 0); /* NOLINT */
                        curpos = curpos.y_rotated(camera_rotation_y).x_rotated(camera_rotation_x); /* is still "at the origin" */
                        line_t ray = line_between(camera_pos, camera_n + curpos); /* t positive is "forward" */
                        wchar_t current_char = L' ';
                        rgb_t current_color;

                        total_bounce_count += scene.render_ray(ray, current_color, current_char, 0.0f, nc);

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
        vthreads.push_back(new std::jthread(render_region, i));
    }

    while (true) {
        while (get_current_time() - last_time < wait_us_per_frame) {;}
        std::uint64_t render_time = get_current_time() - last_time;
        last_time = get_current_time();

        while (!ncplane_mutex.try_lock()) {;}
        ncplane_dim_yx(plane, &dimy, &dimx);
        
        std::uint32_t chin = notcurses_get_nblock(nc, nullptr);
        ncplane_mutex.unlock();
        vec3_t move_d;
        if (chin == L'\t') { break; }
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

        else if (chin == L'j' || chin == L'J') { camera_rotation_x += angle_increment; }
        else if (chin == L'l' || chin == L'L') { camera_rotation_x -= angle_increment; }
        else if ((chin == L'i' || chin == L'I') && camera_rotation_y < std::numbers::pi / 2.0f - angle_increment) {
            camera_rotation_y += angle_increment;
        } else if ((chin == L'k' || chin == L'K') && camera_rotation_y > - (std::numbers::pi / 2.0f - angle_increment)) {
            camera_rotation_y -= angle_increment;
        }

        else if (chin == L' ') {
            gobj_t *closest_obj = nullptr;
            float closest_t = std::numeric_limits<float>::max();
            itimes.clear();
            thisout.clear();
            scene.intersect_ray(line_between(camera_pos, camera_n), nullptr, closest_obj, closest_t, nc);
            if (closest_obj != nullptr) {
                delete closest_obj; /* does not modify closest_obj */
                std::erase_if(scene.objects, [&closest_obj](const gobj_t * const &pgobj) { return pgobj == closest_obj; });
            }
        }

        else if (chin == L't') {
            std::string in;
            while (true) {
                std::uint32_t chin = notcurses_get_blocking(nc, nullptr);
                if (chin == L'\\') {
                    break;
                }
                in.push_back(static_cast<char>(chin));
            }
            base_camera_n.z += std::strtof(in.c_str(), nullptr);
        }

        camera_n = rotated_x_about(rotated_y_about(base_camera_n, camera_pos, camera_rotation_y), camera_pos, camera_rotation_x);

        move_d = move_d.x_rotated(camera_rotation_x);
        camera_pos += move_d;
        camera_n += move_d;
        base_camera_n += move_d;
        std::uint64_t world_objects_count = 0;
        for (const scene_t * const scene : world.scenes) {
            world_objects_count += scene->objects.size();
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
        ncplane_printf_yx(plane, 3, 0, "nx: %-+4.2f ", camera_n.x);
        ncplane_printf_yx(plane, 4, 0, "ny: %-+4.2f ", camera_n.y);
        ncplane_printf_yx(plane, 5, 0, "nz: %-+4.2f ", camera_n.z);
        ncplane_printf_yx(plane, 6, 0, "scene obj: %zu ", scene.objects.size());
        ncplane_printf_yx(plane, 7, 0, "world obj: %zu ", world_objects_count);
        ncplane_printf_yx(plane, 8, 0, "render: %lu µs ", render_time);
        ncplane_printf_yx(plane, 9, 0, "bounces: %lu ", total_bounce_count);
        ncplane_putstr_yx(plane, 10, 0, "x regions: ");
        for (std::int32_t i = 0; i < regions.size(); i++) {
            ncplane_printf_yx(plane, i + 11, 0, "%i : %i ", regions[i].first, regions[i].second);
        }

        notcurses_render(nc);
        ncplane_mutex.unlock();
    }

    /* no idea why stuff gets blocked here, this does not fix */

    /* destructing jthreads and stopping them */
    for (std::jthread *jt : vthreads) {
        jt->request_stop();
        delete jt;
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
