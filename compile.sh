clang++-16 -std=c++20 -DUSE_MULTIMEDIA=none -o rayt src/common.cc src/models.cc src/rayt.cc $(pkg-config --cflags --libs notcurses) -g -O3 -Wall # -fsanitize=address
