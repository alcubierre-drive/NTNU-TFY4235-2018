#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include "TinyPngOut.h"
#include "helpers.h"

#define LMAX 200
// 200
#define MAXNUM 1000
// 1000
#define PSTEP 800
// 200.0
#define LSTEP 10

// generate a matrix of random numbers (squared) that is padded by a leading and 
// a following row of zeros (due to helical boundary conditions)
int ** random_numbers( int N , double p , int SEED ) {

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r,SEED);

    int ** ranreturn = (int **) malloc((N+2)*sizeof( int * ));

    // set the outer values to 0
    ranreturn[0] = (int *) malloc(N*sizeof( int ));
    ranreturn[N+1] = (int *) malloc(N*sizeof( int ));
    for (int j=1; j<N; j++) {
        ranreturn[0][j] = 0;
        ranreturn[N+1][j] = 0;
    }

    // set the inner values to #ran
    for (int i=1; i<N+1; i++) {
        ranreturn[i] = (int *) malloc(N*sizeof( int ));
        for (int j=0; j<N; j++) {
            if (gsl_rng_uniform(r) <= p) {
                ranreturn[i][j] = 1;
            } else {
                ranreturn[i][j] = 0;
            }
        }
    }

    gsl_rng_free(r);
    return ranreturn;
}

// map to use for the matrix:
/*
 * 0: no tree
 * 1: healthy tree
 * 2: lit tree
 * 3: burnt tree
 * 4: ashes
 */

// This implementation uses helical boundary conditions, since the array is
// C-style ordered. The most important function, a time step.
void timestep ( int ** matrix, int N ) {
    // second loop: modify burning status
    for ( int i=1; i<N+1; i++ ) {
        for ( int j=0; j<N; j++ ) {
            int * Pnow = &(matrix[i][j]);
            if (*Pnow == 2) {
                *Pnow = 3;
            } else if (*Pnow == 3) {
                *Pnow = 4;
            }
        }
    }
    // first loop: ignite the trees.
    for ( int i=1; i<N+1; i++ ) {
        for ( int j=0; j<N; j++ ) {
            int left_neighbor[2] = {i+(j-1)/N,(j-1)%N};
            int right_neighbor[2] = {i+(j+1)/N,(j+1)%N};
            int lower_neighbor[2] = {i+1,j};
            if (matrix[i][j] == 1) {
                if (
                        matrix[left_neighbor[0]][left_neighbor[1]] == 2 ||
                        matrix[left_neighbor[0]][left_neighbor[1]] == 3 ||
                        matrix[right_neighbor[0]][right_neighbor[1]] == 2 ||
                        matrix[right_neighbor[0]][right_neighbor[1]] == 3
                    ) {
                    matrix[i][j] = 2;
                }
            } else if (matrix[i][j] == 2 || matrix[i][j] == 3) {
                if (matrix[right_neighbor[0]][right_neighbor[1]] == 1) {
                    matrix[right_neighbor[0]][right_neighbor[1]] = 2;
                }
                if (matrix[lower_neighbor[0]][lower_neighbor[1]] == 1) {
                    matrix[lower_neighbor[0]][lower_neighbor[1]] = 2;
                }
            }
        }
    }
}

// start the fire in the first row.
void start_fire ( int ** matrix, int N) {
    for ( int j=0; j<N; j++ ) {
        if ( matrix[1][j] == 1 ) {
            matrix[1][j] = 2;
        }
    }
}

// function to count types of trees ("which")
int count_stuff ( int ** matrix, int N, int which) {
    int mytrees = 0;
    for (int i=1; i<N+1; i++) {
        for (int j=0; j<N; j++) {
            if (matrix[i][j] == which) {
                mytrees ++;
            }
        }
    }
    return mytrees;
}

// a wrapper that does the whole simulation
void simulation (int size, double p , int SEED, long unsigned int * Presult) {
    int ** my_trees = random_numbers( size, p , SEED );
    int tree_num_start = count_stuff(my_trees, size, 1);

    start_fire(my_trees, size);
    //pprint_matrix(stdout,my_trees,size+2,size);

    int last_val = 0;
    int new_val = 0;
    int t = 2;
    timestep(my_trees, size);
    timestep(my_trees, size);
    do {
        last_val = count_stuff(my_trees, size, 4);
        timestep(my_trees, size);
        t++;
        new_val = count_stuff(my_trees, size, 4);
        if (last_val == new_val) {
            break;
        }
    } while(1);

    free_matrix(my_trees ,size+2);

    *(Presult) += tree_num_start;
    *(Presult+1) += new_val;
    *(Presult+2) += t;
}

// a "wrapper of the wrapper" that saves the images
void save_example_simulation( int size, double p, int SEED , char * fname_prefix ) {

    int ** mymatrix = random_numbers(size, p, SEED);
    char filename[32];
    sprintf(filename,"fire/%s_%03d.png",fname_prefix,0);
    print_image_of_forest( filename, mymatrix, size+2, size);

    start_fire(mymatrix,size);

    sprintf(filename,"fire/%s_%03d.png",fname_prefix,1);
    print_image_of_forest( filename, mymatrix, size+2, size);

    int last_val = 0;
    int new_val = 0;
    int t = 2;
    do {
        last_val = count_stuff(mymatrix, size, 4);
        timestep(mymatrix, size);

        sprintf(filename,"fire/%s_%03d.png",fname_prefix,t);
        print_image_of_forest( filename, mymatrix, size+2, size);

        t++;
        new_val = count_stuff(mymatrix, size, 4);
        if (last_val == new_val) {
            if (count_stuff(mymatrix,size,2) == 0)
                break;
        }
    } while(1);

    for (int k=0; k<10; k++) {
        sprintf(filename,"fire/%s_%03d.png",fname_prefix,t+k);
        print_image_of_forest( filename, mymatrix, size+2, size);
        timestep(mymatrix,size);
    }
    free_matrix(mymatrix, size+2);
}

int main () {
    // here the main program (with really simple parallelization. The result is 
    // just printed to stdout, then piped to a file within the shell.
    printf("#L p 1 4 t\n");
    #pragma omp parallel for
    for (int l = LSTEP; l<=LMAX; l+=LSTEP) {
        for (int i = 1; i<PSTEP; i++) {
            long unsigned int result[3] = {0,0,0};
            double p;
            p = (double)i/(double)PSTEP;
            for (int seed=0; seed<MAXNUM; seed++) {
                simulation( l, p, seed, &(result[0]) );
            }
            printf("%i %.3e %.6e %.6e %.6e\n", l, p,
                    (double)result[0]/(double)MAXNUM,
                    (double)result[1]/(double)MAXNUM,
                    (double)result[2]/(double)MAXNUM);
        }
    }
    /*
    save_example_simulation(900, 0.63, 1, "L900_p0.63");
    save_example_simulation(300, 0.65, 2, "L300_p0.65");
    save_example_simulation(3000, 0.65, 3, "L3000_p0.67");
    */
    // save some images
    save_example_simulation(300, 0.8, 3, "L300_p0.8");
    save_example_simulation(300, 0.4, 3, "L300_p0.4");
    return 0;
}
