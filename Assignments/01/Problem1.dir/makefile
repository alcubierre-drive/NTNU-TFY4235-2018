CC=gcc
CFILE=calcs.c
OFILE=calcs.o
SOFILE=calcs.so
LIBS=-g -lm -lgsl -lgslcblas
#LIBS=-lm -lgsl -lgslcblas -L/usr/X11R6/lib/ -lGL -lGLU -lglut
all:
	$(CC) -O2 -Wall -fPIC -fopenmp -c $(CFILE)
	$(CC) -O2 -Wall -shared $(LIBS) -fopenmp -o $(SOFILE) $(OFILE)
	#@grep -n --color TODO $(CFILE)
