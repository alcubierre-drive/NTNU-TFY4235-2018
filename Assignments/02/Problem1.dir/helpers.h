#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// print a matrix to the terminal, with colors
void pprint_matrix(FILE * file, int ** matrix, int Ni, int Nj) {
    int myval;
    for ( int i=0; i<Ni; i++ ) {
        for ( int j=0; j<Nj; j++) {
            myval = matrix[i][j];
            if (myval == 0) {
                fprintf(file,"%s%s ",ANSI_COLOR_BLUE,"0");
            } else if (myval == 1) {
                fprintf(file,"%s%s ",ANSI_COLOR_GREEN,"1");
            } else if (myval == 2) {
                fprintf(file,"%s%s ",ANSI_COLOR_YELLOW,"2");
            } else if (myval == 3) {
                fprintf(file,"%s%s ",ANSI_COLOR_RED,"3");
            } else if (myval == 4) {
                fprintf(file,"%s%s ",ANSI_COLOR_RESET,"4");
            }
        }
        fprintf(file,"\n");
    }
}

// define some colors for the following function that saves the forest to a .png 
// file using the TinyPngOut library
uint8_t red[3] = {0xFF,0x00,0x00};
uint8_t blue[3] = {0x00,0x00,0xFF};
uint8_t green[3] = {0x00,0xFF,0x00};
uint8_t yellow[3] = {0xFF,0xFF,0x00};
uint8_t grey[3] = {0x77,0x77,0x77};

void print_image_of_forest( char * filename, int ** matrix, int Ni, int Nj) {
    const int height = Ni;
    const int width = Nj;
    uint8_t * pixels = (uint8_t *) malloc(3*Ni*Nj*sizeof(uint8_t));
    for (int i=0; i<Ni; i++) {
        for (int j=0; j<Nj; j++) {
            int myidx = 3*(Nj*i+j);
            int myval = matrix[i][j];
            if (myval == 0) {
                pixels[myidx+0] = blue[0];
                pixels[myidx+1] = blue[1];
                pixels[myidx+2] = blue[2];
            } else if (myval == 1) {
                pixels[myidx+0] = green[0];
                pixels[myidx+1] = green[1];
                pixels[myidx+2] = green[2];
            } else if (myval == 2) {
                pixels[myidx+0] = yellow[0];
                pixels[myidx+1] = yellow[1];
                pixels[myidx+2] = yellow[2];
            } else if (myval == 3) {
                pixels[myidx+0] = red[0];
                pixels[myidx+1] = red[1];
                pixels[myidx+2] = red[2];
            } else if (myval == 4) {
                pixels[myidx+0] = grey[0];
                pixels[myidx+1] = grey[1];
                pixels[myidx+2] = grey[2];
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

// print a matrix
void print_matrix ( FILE * file, int ** matrix, int Ni, int Nj ) {
    for ( int i = 0; i < Ni; i++ ) {
        for (int j = 0; j < Nj; j++ ) {
            fprintf(file,"%i ",matrix[i][j]);
        }
        fprintf(file,"\n");
    }
}

// free the memory of a matrix
void free_matrix ( int ** matrix, int Ni ) {
    for (int i=0; i<Ni; i++) {
        free(matrix[i]);
    }
    free(matrix);
}
