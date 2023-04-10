# rayt
A little project raytracer made from scratch for fun.
Requires notcurses: `sudo apt install libnotcurses-dev` or make and install it from the [repo](https://github.com/dankamongmen/notcurses).
Do not expect any good practices to come out of this code. It may or may not work.


If on your terminal it hangs while starting up, it may be from a race condition in notcurses that happens in some terminals. The solution that has worked for me is to run it using strace like such:
```sh
strace -o /dev/null ./rayt
```
I won't pretend to know why it works, but most of the time it fixes it. Or you could just try another terminal emulator.


Now you can generate your own gradients with `gradient.py`. Put your source chars for the gradient in `gradient_source_chars.txt` as demonstrated, and then run `python gradient.py <font file>` to order the characters by boldness. `<font file>` is a path to the font to use, and you probably want to use whatever font your terminal uses. Works with OpenType (`.otf`) and TrueType (`.ttf`) fonts.
