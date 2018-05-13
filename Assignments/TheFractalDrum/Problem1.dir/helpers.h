#define PI 3.141592653589793

int number_of_points ( int l ) {
    return pow(8,l);
}

void print_limit_memory (int l) {
    int i = sizeof(double);
    double num = 4*pow(8,l);
    int mb = num*i/1024/1024;
    fprintf(stderr, "MINIMAL MEMORY USAGE:\nl = %i: %i MB\n", l, mb);
}

// the grid spacing delta
double delta ( int l ) {
    return 1/pow(8,l);
}

// add two vectors (2D)
void add_vector ( double * pRes, double vec_one[2], double vec_two[2] ) {
    *pRes = vec_one[0]+vec_two[0];
    *(pRes+1) = vec_one[1]+vec_two[1];
}

// substract two vectors (2D)
void diff_vector ( double * pRes, double vec_one[2], double vec_two[2] ) {
    *pRes = vec_one[0]-vec_two[0];
    *(pRes+1) = vec_one[1]-vec_two[1];
}

// multiply a vector (2D) with a scalar
void mult_vector ( double * pRes, double scal ) {
    *pRes *= scal;
    *(pRes+1) *= scal;
}

// rotate a vector by 90 degrees (positive)
void rotate_vector_plus ( double *x, double *y ) {
    double val_x = *x;
    double val_y = *y;
    *x = -val_y;
    *y = val_x;
}

// rotate a vector by 90 degrees (negative)
void rotate_vector_minus ( double *x, double *y ) {
    double val_x = *x;
    double val_y = *y;
    *x = val_y;
    *y = -val_x;
}

// get the angle between two vectors
double angle ( double x1, double x2, double y1, double y2 ) {
    double y = y2-x2;
    double x = y1-x1;
    double ang = atan2(y,x);
    return ang;
}

// this function is not used
double check_for_pi ( double angle_diff ) {
    double result;
    if (angle_diff < -PI) {
        result = angle_diff + 2*PI;
    } else if (angle_diff > PI) {
        result = angle_diff - 2*PI;
    } else {
        result = angle_diff;
    }
    return result;
}
