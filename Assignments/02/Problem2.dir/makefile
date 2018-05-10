CC=gcc
CFILE=calcs.c
EXECFILE=calcs
LIBS=-g -lm -lgsl -lgslcblas -fopenmp
#LIBS=-lm -lgsl -lgslcblas -L/usr/X11R6/lib/ -lGL -lGLU -lglut

all:
	mkdir -p data_plot/ data_time/ img/ plots/
	$(CC) -O2 -Wall $(LIBS) -o $(EXECFILE) TinyPngOut.o $(CFILE)
	#@grep -n --color TODO $(CFILE)

# just to do the vim-stuff
types: types.vim
types.vim: *.[ch]
	ctags --c-kinds=gstu -o- *.[ch] |\
		awk 'BEGIN{printf("syntax keyword Type\t")}\
		{printf("%s ", $$1)}END{print ""}' > $@
