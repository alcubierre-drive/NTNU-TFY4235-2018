CFILE=calcs.c
OFILE=calcs.o
SOFILE=calcs.so

all:
	gcc -O2 -Wall -Wextra -pedantic -fPIC -fopenmp -c $(CFILE)
	gcc -O2 -shared -lgsl -lgslcblas -lm -fopenmp -o $(SOFILE) $(OFILE)
