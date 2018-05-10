#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include "helpers.h"

const int LMAX = 4;
const int OFFGRID = 1;
const double ANGLE_TOLERANCE = 0.01;
// REAL MAXIMUM ON THIS MACHINE IS 8!!!

// get new points of the fractal, used "recursively" later on.
void gen_new_points ( double * pRes, double vec_one[2], double vec_two[2] ) {
    double * new_points = (double *) malloc(2*sizeof(double));
    // calculate distance vector
    double * dist = (double *) malloc(2*sizeof(double));
    diff_vector( &dist[0], vec_two, vec_one );
    mult_vector( &dist[0], 0.25);
    // do the steps of the fractal
    // first vector is the same
    *pRes = vec_one[0];
    *(pRes+1) = vec_one[1];
    // second vector
    add_vector ( pRes+2, vec_one, dist );
    // third vector
    rotate_vector_plus ( &dist[0], &dist[1] );
    new_points[0] = *(pRes+2);
    new_points[1] = *(pRes+3);
    add_vector ( pRes+4, new_points, dist );
    // fourth vector
    rotate_vector_minus ( &dist[0], &dist[1] );
    new_points[0] = *(pRes+4);
    new_points[1] = *(pRes+5);
    add_vector ( pRes+6, new_points, dist );
    // fifth vector
    rotate_vector_minus ( &dist[0], &dist[1] );
    new_points[0] = *(pRes+6);
    new_points[1] = *(pRes+7);
    add_vector ( pRes+8, dist, new_points );
    // sixth vector
    new_points[0] = *(pRes+8);
    new_points[1] = *(pRes+9);
    add_vector ( pRes+10, dist, new_points );
    // seventh vector
    rotate_vector_plus ( &dist[0], &dist[1] );
    new_points[0] = *(pRes+10);
    new_points[1] = *(pRes+11);
    add_vector ( pRes+12, dist, new_points );
    // eighth vector
    rotate_vector_plus ( &dist[0], &dist[1] );
    new_points[0] = *(pRes+12);
    new_points[1] = *(pRes+13);
    add_vector ( pRes+14, dist, new_points );
    free(dist);
}

// function to generate the vector of one side of the boundary
double ** genvec ( const int num, const int l) {
    double ** returnarray = (double **) malloc((num+1)*sizeof( double * ));
    for ( int i = 0; i < num+1; i ++ ) {
        returnarray[i] = (double *) malloc(2*sizeof( double ));
    }
    returnarray[num][0] = 1.0;
    for ( int i = 0; i < l; i++ ) {
        int idx_delta = num/pow(8,i);
        for ( int j = 0; j < pow(8,i); j++ ) {
            double * tmpvec = (double *) malloc(16*sizeof(double));
            gen_new_points(
                    &tmpvec[0],
                    returnarray[idx_delta*j],
                    returnarray[idx_delta*(j+1)]
            );
            for ( int k = 0; k < 8; k++ ) {
                returnarray[idx_delta*j+idx_delta/8*k][0] = tmpvec[2*k];
                returnarray[idx_delta*j+idx_delta/8*k][1] = tmpvec[2*k+1];
            }
        }
    }
    return returnarray;
}

// helper function to transform one side of the boundary into another
double ** shift_big_vec ( double ** big_vec, const int num, double shift[2] ) {
    double ** returnarray = (double **) malloc(num*sizeof( double * ));
    for ( int i = 0; i < num; i ++ ) {
        returnarray[i] = (double *) malloc(2*sizeof( double ));
        returnarray[i][0] = big_vec[i][0] + shift[0];
        returnarray[i][1] = big_vec[i][1] + shift[1];
    }
    return returnarray;
}

// helper function to transform one side of the boundary into another
double ** rotate_big_vec ( double ** big_vec, const int num, int dir ) {
    double ** returnarray = (double **) malloc(num*sizeof( double * ));
    if (dir == 0) {
        // negative
        for ( int i = 0; i < num; i++ ) {
            returnarray[i] = (double *) malloc(2*sizeof(double));
            returnarray[i][0] = big_vec[i][1];
            returnarray[i][1] = -big_vec[i][0];
        }
    } else {
        // positive
        for ( int i = 0; i < num; i++ ) {
            returnarray[i] = (double *) malloc(2*sizeof(double));
            returnarray[i][0] = -big_vec[i][1];
            returnarray[i][1] = big_vec[i][0];
        }
    }
    return returnarray;
}

// final function to create the whole boundary for a given generation l
double ** create_boundary ( const int l ) {
    // the final array's length will be 4*num.
    int num = number_of_points(l);
    double ** base = genvec(num,l);
    double shift[2] = {0.0,1.0};
    double ** upper = shift_big_vec( base, num, shift );
    double ** tmpright = rotate_big_vec( upper, num, 0 );
    double ** right = shift_big_vec( tmpright, num, shift );
    double ** left = rotate_big_vec( base, num, 1 );
    double ** all_boundaries = (double **) malloc(2*sizeof( double * ));
    all_boundaries[0] = (double *) malloc(4*num*sizeof(double));
    all_boundaries[1] = (double *) malloc(4*num*sizeof(double));
    for ( int i = 0; i < num; i++ ) {
        //all_boundaries[i] = (double *) malloc(2*sizeof(double));
        all_boundaries[0][i] = base[num-i][0];
        all_boundaries[1][i] = base[num-i][1];
        //free(base[num-i]);
    }
    free(base);
    for ( int i = num; i < 2*num; i++ ) {
        //all_boundaries[i] = (double *) malloc(2*sizeof(double));
        all_boundaries[0][i] = left[i-num][0];
        all_boundaries[1][i] = left[i-num][1];
        //free(left[i-num]);
    }
    free(left);
    for ( int i = 2*num; i < 3*num; i++ ) {
        //all_boundaries[i] = (double *) malloc(2*sizeof(double));
        all_boundaries[0][i] = upper[i-2*num][0];
        all_boundaries[1][i] = upper[i-2*num][1];
        //free(upper[i-2*num]);
    }
    free(upper);
    for ( int i = 3*num; i < 4*num; i++ ) {
        //all_boundaries[i] = (double *) malloc(2*sizeof(double));
        all_boundaries[0][i] = right[i-3*num][0];
        all_boundaries[1][i] = right[i-3*num][1];
        //free(right[i-3*num]);
    }
    free(right);
    return all_boundaries;
}

// not used
int index_converter( int i, int j, const int GRIDSIZE ) {
    int flat_idx = i*GRIDSIZE + j;
    return flat_idx;
}

// not used
int map_grid_to_final_vec( int i, int j, const int GRIDSIZE,
        int * IN_MATRIX ) {
    int flat_idx = index_converter(i,j,GRIDSIZE);
    int final_idx = 0;
    for ( int k = 0; k < flat_idx; k++ ) {
        final_idx += IN_MATRIX[k];
    }
    return final_idx - 1;
}
