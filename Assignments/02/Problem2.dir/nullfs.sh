#!/bin/bash

mkdir DEVNULL/
nullfs DEVNULL/
export OMP_NUM_THREADS=8
./calcs
fusermount -u DEVNULL/
rmdir DEVNULL/

# some ugly "hack" to not save all the data when generating images of the
# simulation. Uses a package "nullfs" that needs to be installed on the system.
# The Arch Linux AUR snapshot can be found in Libraries/.
