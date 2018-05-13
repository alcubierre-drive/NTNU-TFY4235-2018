SDL := $(shell pkg-config --cflags --libs sdl2)
SDL := -lSDL2_mixer $(SDL)
GLUT=-lGL -lGLU -lglut
LIBNOTIFY := $(shell pkg-config --cflags --libs libnotify)
OTHER=-fopenmp -lm
CFILE=reaction.c

all:
	mkdir -p tests/ const_sim/ var_sim/ st_sim/ var_sim_schemes/ D_sim/
	$(CC) -O3 -o reaction $(CFILE) $(OTHER)
