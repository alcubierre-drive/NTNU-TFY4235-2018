CC=gcc
CFILE=calcs.c
EXECFILE=calcs
LIBFILE=TinyPngOut.o
LIBS=-g -lm -lgsl -lgslcblas -fopenmp
#LIBS=-lm -lgsl -lgslcblas -L/usr/X11R6/lib/ -lGL -lGLU -lglut
all:
	$(CC) -O2 -Wall $(LIBS) -o $(EXECFILE) $(LIBFILE) $(CFILE)
	mkdir -p plots/
	@echo "set the environment variable OMP_NUM_THREADS to e.g. 7"
	@#grep -n --color TODO $(CFILE)
