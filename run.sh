# do this to force notcurses to actually use rgb escape codes
# if this does not work and ends up with weird output, just comment it out
export COLORTERM=truecolor

sudo false # to make sure we have sudo perms

(sleep 2 && sudo killall -s SIGKILL strace 2>/dev/null) &
strace -f -o /dev/null -D ./rayt -f480 -j12 -s2000
