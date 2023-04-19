# do this to force notcurses to actually use rgb escape codes
# if this does not work and ends up with weird output, just comment it out
export COLORTERM=truecolor

(sleep 2 && sudo killall -s SIGKILL strace) &
strace -f -o /dev/null -D ./rayt
pidwait rayt
