clang++ -o rayt src/rayt.cc $(pkg-config --cflags --libs notcurses) -g -std=c++20 -Wall
