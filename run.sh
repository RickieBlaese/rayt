# do this to force notcurses to actually use rgb escape codes
# if this does not work and ends up with weird output, just comment it out
export COLORTERM=24bit

strace -o /dev/null ./rayt
