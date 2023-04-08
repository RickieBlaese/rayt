# rayt
A little project raytracer made from scratch for fun.
Requires notcurses: `sudo apt install libnotcurses-dev` or make and install it from the [repo](https://github.com/dankamongmen/notcurses).
Do not expect any good practices to come out of this code. It may or may not work.


If on your terminal it hangs while starting up, it may be from a race condition in notcurses that happens in some terminals. The solution that has worked for me is to run it using strace like such:
```sh
strace -o /dev/null ./rayt
```
I won't pretend to know why it works, but most of the time it fixes it. Or you could just try another terminal emulator.
