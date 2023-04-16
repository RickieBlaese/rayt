clang++ -o rayt src/common.cc src/models.cc src/rayt.cc $(pkg-config --cflags --libs notcurses) -g -std=c++20 -Wall -O3
