#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include "TinyPngOut.h"
#include "helpers.h"

// generate a random 0-1 int-square
int ** random_numbers_square( int N , int SEED ) {

    // the infected 1% is hardcoded. reduces errors.

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r,SEED);

    int ** ranreturn = (int **) malloc(N*sizeof( int * ));

    // set the inner values to #ran
    for (int i=0; i<N; i++) {
        ranreturn[i] = (int *) malloc(N*sizeof( int ));
        for (int j=0; j<N; j++) {
            if (gsl_rng_uniform(r) <= 0.01) {
                ranreturn[i][j] = 1;
            } else {
                ranreturn[i][j] = 0;
            }
        }
    }

    gsl_rng_free(r);
    return ranreturn;
}

// generate the infected population
individual ** infect_population( int N, int SEED ) {
    int ** infection = random_numbers_square( N, SEED );
    individual ** result = (individual **) malloc(N*sizeof(individual * ));
    for (int i=0; i<N; i++) {
        result[i] = (individual *) malloc(N*sizeof(individual));
        for (int j=0; j<N; j++) {
            result[i][j].sick_current = infection[i][j];
            result[i][j].sick_old = 0;
            result[i][j].pathogen = 0;
            result[i][j].recovered[0] = 0;
            for (int k=1; k<NUMBER_OF_MUTATIONS; k++) {
                result[i][j].recovered[k] = 0;
            }
        }
    }
    free_matrix(infection,N);
    return result;
}

// transmit the current pathogen to the nearest neighbors
void transmit_to_neighbors( individual ** inds, index * neighbors,
        int mutation ) {
    for (int l=0; l<4; l++) {
        inds[neighbors[l].i][neighbors[l].j].pathogen = mutation;
    }
}

// one timestep of the simulation
void timestep ( individual ** inds, int N, double p, double lambda, double q,
        gsl_rng * r) {
    // make current sickness old, and transmit the pathogen to the neighbors
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            index idx; idx.i = i; idx.j = j;
            index * neighbors = malloc(4*sizeof(index));
            // transmit the pathogen
            nearest_neighbors( idx, N , neighbors);
            if (inds[i][j].sick_current != 0) {
                transmit_to_neighbors(inds,neighbors,inds[i][j].sick_current);
            }
            free(neighbors);
            // make current sickness old
            inds[i][j].sick_old = inds[i][j].sick_current;
            inds[i][j].sick_current = 0;
        }
    }

    // infect the nodes, and mutate the pathogen.
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (inds[i][j].pathogen != 0) {
                // do the mutations.
                if (lambda != 0.0) {
                    double lambda_rand = gsl_rng_uniform(r);
                    if (lambda_rand <= lambda) {
                        inds[i][j].pathogen = random_from_double( lambda_rand, lambda );
                        // This is "undefined behaviour":
                        // inds[i][j].pathogen = double_as_llint( &lambda_rand );
                    }
                }
                // do the infections.
                if (inds[i][j].recovered[inds[i][j].pathogen-1] == 1) {
                    // if it has recovered, infect with p*q
                    if (gsl_rng_uniform(r) <= p*q) {
                        inds[i][j].sick_current = inds[i][j].pathogen;
                    }
                } else {
                    // if it has never had that mutation before, p!
                    if (gsl_rng_uniform(r) <= p) {
                        inds[i][j].sick_current = inds[i][j].pathogen;
                    }
                }
                // remove pathogen after this step.
                inds[i][j].pathogen = 0;
            }
        }
    }

    // recover individuals that have not been reinfected.
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (inds[i][j].sick_old != 0 && inds[i][j].sick_current == 0) {
                inds[i][j].recovered[inds[i][j].sick_old-1] = 1;
            }
        }
    }
}

// count the current infections
int count_infections (individual ** inds, int N) {
    int result = 0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (inds[i][j].sick_current != 0) {
                result++;
            }
        }
    }
    return result;
}

// count the one-time infections
int count_one_time_infections(individual ** inds, int N) {
    int result = 0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            for (int k=0; k<NUMBER_OF_MUTATIONS; k++) {
                if (inds[i][j].sick_current != 0 || inds[i][j].recovered[k] != 0) {
                    result++;
                    break;
                }
            }
        }
    }
    return result;
}

// generate the densities, and eventually write some images
double ** rho ( individual ** inds, int N, double Nsquare, double p,
        double lambda, double q, gsl_rng * r) {
    int not_too_many;
    if (MAXTIME < 500) {
        not_too_many = 1;
    } else {
        not_too_many = 0;
    }
    double ** result = malloc(2*sizeof(double*));
    result[0] = malloc(MAXTIME*sizeof(double));
    result[1] = malloc(MAXTIME*sizeof(double));
    for (int t=0; t<MAXTIME; t++) {
        result[0][t] = (double)count_infections( inds, N )/Nsquare;
        result[1][t] = (double)count_one_time_infections( inds, N )/Nsquare;
        timestep( inds, N, p, lambda, q, r);
        if (PRINT_IMAGES == 1 && not_too_many == 1) {
            char filename[128];
            sprintf(filename,"%sp%.3fq_%.3f_lambda%.3f_t%05d.png","img/",p,q,lambda,t+1);
            print_image_of_illness(filename, inds, N, N);
        }
    }
    return result;
}

// generate the densities for the case that the extinction time should be 
// determined -- breaks the loop after it is certain that it's done. Also writes 
// some images
double ** rho_extinction (individual ** inds, int N, double Nsquare, double p,
        gsl_rng * r, int * tmax) {
    if (!(p < 1.)) {
        exit(1);
    } else {
        double ** result = malloc(2*sizeof(double*));
        result[0] = malloc(MAXTIME*sizeof(double));
        result[1] = malloc(MAXTIME*sizeof(double));
        for (int t=0; t<MAXTIME; t++) {
            result[0][t] = (double)count_infections( inds, N )/Nsquare;
            result[1][t] = (double)count_one_time_infections( inds, N )/Nsquare;
            if (result[0][t] < 0.5/(double)Nsquare) {
                *tmax = t;
                break;
            }
            timestep( inds, N, p, 0.0, 0.0, r);
        }
        for (int t=*tmax;t<MAXTIME;t++) {
            result[0][t] = 0.0;
            result[1][t] = result[1][*tmax];
        }
        return result;
    }
}

// basic function to get *one* simulation done. Can either be with the break 
// after the extinction or the regular one.
void simulation ( char * fname, double p, double q, double lambda,
        gsl_rng * r, int SEED, int MODE ) {
        if (MODE == 1 && !(p<1.)) {
            exit(1);
        }
        FILE * myfile = fopen(fname,"w");
        individual ** inds = infect_population( LATTICE_SIZE, SEED );
        double ** densities;
        int tmax = 0;
        if (MODE == 0) {
            densities = rho(inds, LATTICE_SIZE, pow(LATTICE_SIZE,2), p,
                    lambda, q, r);
        } else if (MODE == 1) {
            densities = rho_extinction(inds, LATTICE_SIZE, pow(LATTICE_SIZE,2),
                    p, r, &tmax);
        } else {
            fprintf(stderr,"No such ´MODE´, use ´1´ or ´2´\n");
            exit(1);
        }
        if (MODE == 1) {
            fprintf(myfile,"# tmax = %i\n", tmax);
        }
        for (int k=0; k<MAXTIME; k++) {
            fprintf(myfile,"%i %.7f %.7f\n",k,densities[0][k],densities[1][k]);
        }

        free(densities[0]); free(densities[1]); free_individuals(inds,LATTICE_SIZE);
        fclose(myfile);
}

// wrapper to generate (and write) some data of the simulatins in ranges of 
// p_range and q_range.
void generate_data_plots ( char * prefix, range p_range, range q_range,
        double lambda, gsl_rng * r, int SEED) {
#pragma omp parallel for
    for (int i=0; i<p_range.N; i++) {
        double my_p = range_i(p_range,i);
        for (int j=0; j<q_range.N; j++) {
            double my_q = range_i(q_range,j);
            char filename[64];
            sprintf(filename,"%sp%.3f_q%.3f_lambda%.3f.dat",prefix,my_p,my_q,lambda);
            // MODE = 0
            simulation( filename, my_p, my_q, lambda, r, SEED, 0 );
        }
    }
}

// wrapper to generate some data for extinction_time mode, for p in p_range, 
// takes the mean over all seeds specified in seedrange. Only writes t_max and 
// rho_2.
void generate_time_data ( char * prefix, range p_range,
        gsl_rng * r, range seedrange) {
    double * rho2 = malloc(p_range.N*sizeof(double));
    double * tmax = malloc(p_range.N*sizeof(double));
#pragma omp parallel for
    for (int i=0; i<p_range.N; i++) {
        double my_p = range_i(p_range,i);
        // averages
        rho2[i] = 0.0;
        tmax[i] = 0.0;
        for (int seed = (int)seedrange.min; seed <= (int)seedrange.max; seed++) {
            individual ** inds = infect_population(LATTICE_SIZE, seed);
            int currtmax;
            double ** currdens = rho_extinction(inds, LATTICE_SIZE,
                    pow(LATTICE_SIZE,2), my_p, r, &currtmax);
            rho2[i] += currdens[1][currtmax]/(double)seedrange.N;
            tmax[i] += currtmax/(double)seedrange.N;
            free(currdens[0]); free(currdens[1]); free_individuals(inds,LATTICE_SIZE);
        }
    }
    // writeout
    char filename[64];
    sprintf(filename,"%sp%.3f-%.3f-%05d.dat",prefix,p_range.min,p_range.max,p_range.N);
    FILE * myfile = fopen(filename,"w");
    for (int i=0; i<p_range.N; i++) {
        double my_p = range_i(p_range,i);
        fprintf(myfile,"%.7f %.7f %.7f\n",my_p,rho2[i],tmax[i]);
    }
    free(rho2); free(tmax); fclose(myfile);
}

int main () {
    // preparation
    int initial_seed = 1;
    int other_seed = 2;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r,other_seed);

    // this time too lazy to write some logic to handle this. So for each 
    // dataset produced, the code needs to be changed accordingly (just comment 
    // / uncomment and change constants defined *in helpers.h* or parameters 
    // defined *directly* here)
    range p, q;
    p.min = 0.1; p.max = 1.01; p.N = 100;
    q.min = 0.0; q.max = 1.5; q.N = 3;
    for (double lambda = 0.0; !(lambda > 1.0); lambda += 0.5) {
        // UNCOMMENT IF WANT TO GENERATE
        // generate_data_plots("data_plot/", p, q, lambda, r, initial_seed);
    }

    p.min = 0.0; p.max = 1.0; p.N = 999;
    range sran;
    sran.min = 0; sran.max = 1000; sran.N = 1000;
    // UNCOMMENT IF WANT TO GENERATE
    // generate_time_data("data_time/", p, r, sran);

    // The following snippet will just save images for the configurations in the 
    // array configs[][]. Thus, the prefix which is given to the simulation() 
    // function is the folder ´DEVNULL/´ which has a fuse-filesystem mounted 
    // that points to /dev/null (fast discard of the data) -- this is done 
    // because it was the fastest solution; and works.
    // The Arch Linux PKGBUILD (from the AUR) to install the package nullfs can 
    // be found in the gzipped tar archive nullfs.tar.gz.
    double configs[8][3] = {{0.5,0.5,0.0},{0.5,0.0,0.0},{1.0,0.5,0.0},{1.0,0.0,0.0},
        {0.5,0.5,0.2},{0.5,0.0,0.2},{1.0,0.5,0.2},{1.0,0.0,0.2}};
#pragma omp parallel for
    for (int i=0; i<8; i++) {
        simulation ( "DEVNULL/",
                configs[i][0], configs[i][1], configs[i][2],
                r, initial_seed, 0 );
    }

    gsl_rng_free(r);
    return 0;
}
