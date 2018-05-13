#define NUMBER_OF_MUTATIONS (int)100
#define MAXTIME (int)200
#define LATTICE_SIZE (int)500

#define PRINT_IMAGES (int)1

// define the datatype of an index
struct index {
    int i;
    int j;
};
typedef struct index index;

// define the datatype of an individual. It carries the information about what 
// pathogen it has, if it was sick last timestep, if it is sick this timestep 
// and for which of the mutations it has recovered.
struct individual {
    int pathogen;
    int sick_old;
    int sick_current;
    int recovered[NUMBER_OF_MUTATIONS];
};
typedef struct individual individual;

// define the datatype of a range, its minimum, maximum and the number of points 
// that should be used.
struct range {
    double min;
    double max;
    int N;
};
typedef struct range range;

// helper to get the current value of a range-type
double range_i( range r, int i ) {
    double result = r.min + (r.max-r.min)/(double)r.N*(double)i;
    return result;
}

// function to get the random part of a double "lambda_rand" smaller than the 
// maximum given "lambda_max". Used in order to save time
int random_from_double( double lambda_rand, double lambda_max ) {
    int result = (int)(lambda_rand/lambda_max*NUMBER_OF_MUTATIONS);
    return result % NUMBER_OF_MUTATIONS + 1;
}

/* A nice "hack" to use less computational resources. But it is not
 * reporducible ("undefined behaviour"), so it won't be used.
long long int double_as_llint( double * d ) {
    long long int * result;
    result = d;
    return (*result) % NUMBER_OF_MUTATIONS + 1;
}
*/

// define some colors to generate the images
uint8_t red[3] = {0xFF,0x00,0x00};
uint8_t blue[3] = {0x00,0x00,0xFF};
uint8_t green[3] = {0x00,0xFF,0x00};
uint8_t yellow[3] = {0xFF,0xFF,0x00};
uint8_t grey[3] = {0x77,0x77,0x77};

// helper to calculate the color that lies 0<x<1 in between of two other colors
uint8_t mean ( uint8_t a, uint8_t b, double f ) {
    uint8_t result = a+(uint8_t)((double)(b-a)*f);
    return result;
}

// get the color array of the "mean" color
void color_transition ( uint8_t * color_one, uint8_t * color_two,
        double between, uint8_t * result ) {
    for (int i=0; i<3; i++) {
        result[i] = mean(color_one[i], color_two[i], between);
    }
}

// generate an image of the current state
void print_image_of_illness( char * filename, individual ** matrix, int Ni, int Nj) {
    const int height = Ni;
    const int width = Nj;
    uint8_t * pixels = (uint8_t *) malloc(3*Ni*Nj*sizeof(uint8_t));
    for (int i=0; i<Ni; i++) {
        for (int j=0; j<Nj; j++) {
            int myidx = 3*(Nj*i+j);
            individual myval = matrix[i][j];
            if (myval.sick_current == 0) {
                pixels[myidx+0] = blue[0];
                pixels[myidx+1] = blue[1];
                pixels[myidx+2] = blue[2];
            } else {
                uint8_t color[3];
                double between = (double)(myval.sick_current-1)/(double)NUMBER_OF_MUTATIONS;
                color_transition(red, yellow, between, color);
                pixels[myidx+0] = color[0];
                pixels[myidx+1] = color[1];
                pixels[myidx+2] = color[2];
            }
        }
    }
    FILE * fout = fopen(filename,"wb");
    struct TinyPngOut pngout;
    TinyPngOut_init(&pngout, fout, width, height);
    TinyPngOut_write(&pngout, pixels, width * height);
    TinyPngOut_write(&pngout, NULL, 0);
    fclose(fout);
    free(pixels);
}

// free a matrix of ints
void free_matrix ( int ** matrix, int Ni ) {
    for (int i=0; i<Ni; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// free a matrix of individuals
void free_individuals ( individual ** inds, int Ni ) {
    for (int i=0; i<Ni; i++) {
        free(inds[i]);
    }
    free(inds);
}

// get the periodical-boundary-conditions-version of an index of a matrix
index periodical_index ( index non_periodical_index , int N ) {
    index result;
    result.i = ((non_periodical_index.i % N) + N) % N;
    result.j = ((non_periodical_index.j % N) + N) % N;
    return result;
}

// print the nearest neighbors of an index to an array of indices
void nearest_neighbors ( index current_index, int N, index * writeout) {
    index run_index;
    run_index = current_index;
    run_index.i += 1;
    writeout[0] = periodical_index( run_index, N );
    run_index.i -= 2;
    writeout[1] = periodical_index( run_index, N );
    run_index.i += 1;
    run_index.j += 1;
    writeout[2] = periodical_index( run_index, N );
    run_index.j -= 2;
    writeout[3] = periodical_index( run_index, N );
}
