#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const int MAXTIME = 100000;
const int MAXDATA = 300000;

double gauss_rannum ( gsl_rng * r ) {
    return gsl_ran_gaussian(r,1.0);
}

double uniform_rannum ( gsl_rng * r ) {
    return -1.0+2.0*gsl_rng_uniform(r);
}

int * get_times_array ( char const type ) {

    if (MAXTIME >= INT_MAX || MAXDATA >= INT_MAX) {
        fprintf(stderr,"ERROR: ´MAXTIME´ or ´MAXDATA´ too big.");
        exit(42);
    }

    int * times_vec = (int *) malloc(MAXDATA*sizeof( int ));

    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    double (*randomP)(gsl_rng *);
    switch (type) {
        case 'g':
            randomP = &gauss_rannum;
            break;
        case 'u':
            randomP = &uniform_rannum;
            break;
        default:
            randomP = &uniform_rannum;
            break;
    }

    int thread_id, nloops;
    #pragma omp parallel private(thread_id, nloops)
    {
        nloops = 0;
        #pragma omp for
        for (int i = 0; i<MAXDATA; i++) {
            double tmpsum = 0;
            int t;
            for (t = 0; t<MAXTIME; t++) {
                double rannum;
                rannum = (*randomP)(r);
                if ( t == 0 ) {
                    tmpsum += fabs(rannum);
                }
                else {
                    tmpsum += rannum;
                    if ( tmpsum <= 0.0 ) {
                        times_vec[i] = t;
                        break;
                    }
                }
            }
            if ( t == MAXTIME-1 ) {
                times_vec[i] = 0;
            }
            ++nloops;
        }
        thread_id = omp_get_thread_num();
        fprintf(stderr, "Thread %d: %d iterations.\n",thread_id, nloops);
    }

    gsl_rng_free(r);
    return times_vec;
}
