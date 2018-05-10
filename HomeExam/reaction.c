#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>

#include "schemes.c"
// Here, the additional testing schemes for variable D are included. They 
// require some functions defined below and thus contain some prototype 
// expressions.

int shift_idx( int idx, int shift, int N ) {
    // shifts an index with shift and then checks bounds
    int result = (idx + N + shift) % N;
    return result;
}

double heaviside( double x ) {
    // just a definition of the heaviside function. returns double.
    if (x >= 0.) {
        return 1.;
    } else {
        return 0.;
    }
}

void save_array_to_file( double ** final, double * xvec, char * fname, int N ) {
    // saves an array of given shape (plus xvec) to a simple text file
    FILE * myfile = fopen(fname, "w");
    fprintf(myfile,"# x, a, b, c, s\n");
    for (int i=0; i<N; i++) {
        fprintf(myfile,"%.6e %.6e %.6e %.6e %.6e\n",
                xvec[i], final[0][i], final[1][i], final[2][i], final[3][i]);
    }
    fclose(myfile);
}

void coupled_constant_single( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    /*
     * params = {da, db, dc, r, n1, n2, c0, b0/a0}
     * in / out = {alpha, beta, gamma, sigma}
     * One time step for the reaction equations. Takes the vectors alpha, beta,
     * gamma and sigma and performs the scheme.
     */

    // some calculations for the parameters.
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    // the actual scheme applied to all four vectors.
    for (int i=0; i<N; i++) {
        out[0][i] = in[0][i] + d_a*( in[0][shift_idx(i,1,N)] +
                in[0][shift_idx(i,-1,N)] - 2.*in[0][i] ) - r*in[0][i]*in[1][i];
        out[1][i] = in[1][i] + d_b*( in[1][shift_idx(i,1,N)] +
                in[1][shift_idx(i,-1,N)] - 2.*in[1][i] ) - r*in[0][i]*in[1][i];
        out[2][i] = in[2][i] + d_c*( in[2][shift_idx(i,1,N)] +
                in[2][shift_idx(i,-1,N)] - 2.*in[2][i] ) + r*in[0][i]*in[1][i] -
                n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }

    // the boundary conditions.
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

double D_func ( double s, double s0 ) {
    // small helper for the diffusivity functional
    double result = 1./(1.+s/s0);
    return result;
}

double D_func_d ( double s, double s0 ) {
    // small helper for the prefactor of the derivative of the diffusivity 
    // functional
    double result = -1./pow(1.+s/s0,2)/s0;
    return result;
}

void coupled_variable_single ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    /*
     * params = {da, db, dc, r, n1, n2, c0, b0/a0, sigma0}
     * in / out = {alpha, beta, gamma, sigma}
     * One time step for the reaction equations. Takes the vectors alpha, beta,
     * gamma and sigma and performs the scheme.
     */

    // some calculations for the parameters.
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    // specific part of the non-constant diffusivity (NOT the same as s0)
    double sigma_0 = params[8];

    // the actual scheme (D_g not constant!) applied to all four vectors.
    for (int i=0; i<N; i++) {
        // each step only one value for D_g
        double D_g = D_func( in[3][i], sigma_0 );
        // and one value for the factor in front of the derivative
        double D_g_diff = D_func_d( in[3][i], sigma_0 ) / 2.;

        out[0][i] = in[0][i] + D_g * d_a*( in[0][shift_idx(i,1,N)] +
                in[0][shift_idx(i,-1,N)] - 2.*in[0][i] ) - r*in[0][i]*in[1][i] +
                D_g_diff * d_a*( (in[3][i]-in[3][shift_idx(i,-1,N)]) * 
                        (in[0][shift_idx(i,1,N)] - in[0][i]) + 
                        (in[3][shift_idx(i,1,N)] - in[3][i]) * 
                        (in[0][i] - in[0][shift_idx(i,-1,N)]) );
        out[1][i] = in[1][i] + D_g * d_b*( in[1][shift_idx(i,1,N)] +
                in[1][shift_idx(i,-1,N)] - 2.*in[1][i] ) - r*in[0][i]*in[1][i] +
                D_g_diff * d_a*( (in[3][i]-in[3][shift_idx(i,-1,N)]) * 
                        (in[1][shift_idx(i,1,N)] - in[1][i]) + 
                        (in[3][shift_idx(i,1,N)] - in[3][i]) * 
                        (in[1][i] - in[1][shift_idx(i,-1,N)]) );
        out[2][i] = in[2][i] + D_g * d_c*( in[2][shift_idx(i,1,N)] +
                in[2][shift_idx(i,-1,N)] - 2.*in[2][i] ) + r*in[0][i]*in[1][i] -
                n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i] +
                D_g_diff * d_a*( (in[3][i]-in[3][shift_idx(i,-1,N)]) * 
                        (in[2][shift_idx(i,1,N)] - in[2][i]) + 
                        (in[3][shift_idx(i,1,N)] - in[3][i]) * 
                        (in[2][i] - in[2][shift_idx(i,-1,N)]) );
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }

    // the boundary conditions.
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

void coupled_stoichio_single( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N, double (*rfun)(double, double) ) {
    /*
     * params = {da, db, dc, r, n1, n2, c0, b0/a0, stoichio_a, stoichio_b}
     * in / out = {alpha, beta, gamma, sigma}
     * One time step for the reaction equations (with changed stoichiometry, 
     * task h)).
     * Takes the vectors alpha, beta, gamma and sigma and performs the scheme.
     *
     * The function rfun is used to change the appearance of R(a,b). it gets 
     * multiplied by r, so it has to be "normalized".
     */

    // some calculations for the parameters.
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    // the stoichiometric coefficient
    double stoichio_a = params[8];
    double stoichio_b = params[9];

    // the actual scheme applied to all four vectors.
    for (int i=0; i<N; i++) {
        out[0][i] = in[0][i] + d_a*( in[0][shift_idx(i,1,N)] +
                in[0][shift_idx(i,-1,N)] - 2.*in[0][i] ) -
                stoichio_a * r * (*rfun)(in[0][i],in[1][i]);
        out[1][i] = in[1][i] + d_b*( in[1][shift_idx(i,1,N)] +
                in[1][shift_idx(i,-1,N)] - 2.*in[1][i] ) -
                stoichio_b * r * (*rfun)(in[0][i],in[1][i]);
        out[2][i] = in[2][i] + d_c*( in[2][shift_idx(i,1,N)] +
                in[2][shift_idx(i,-1,N)] - 2.*in[2][i] ) +
                r * (*rfun)(in[0][i],in[1][i]) -
                n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }

    // the boundary conditions.
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

void coupled_constant( double * params, int N_chi, double d_tau, int N_tau,
        double ** initial, char * fname ) {
    /*
     * params the same as in coupled_constant_single(), only an initial value 
     * for the configuration is needed. This is done to make several images (in 
     * time) of only one calculation, i.e. continue the calculation after saving 
     * the data.
     */

    // calculate d_chi
    double d_chi = 1./((double)N_chi-1.);

    // some memory stuff, also create xvec.
    double * xvec = malloc(N_chi*sizeof(double));
    double ** final = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        final[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        xvec[i] = (double)i/((double)N_chi-1);
    }

    // solve the pde N_tau times, each step: use the result as the initial 
    // condition in the next step.
    for (int t=0; t<N_tau; t++) {
        coupled_constant_single( params, d_tau, d_chi, initial, final, N_chi );
        for (int i=0; i<N_chi; i++) {
            for (int j=0; j<4; j++) {
                initial[j][i] = final[j][i];
            }
        }
    }

    // save data to file.
    save_array_to_file(initial, xvec, fname, N_chi);

    // free the vector final.
    for (int j=0; j<4; j++) {
        free(final[j]);
    }
    free(final);
}

void coupled_variable( double * params, int N_chi, double d_tau, int N_tau,
        double ** initial, char * fname, int MODE) {
    /*
     * params the same as in coupled_variable_single(), only an initial value 
     * for the configuration is needed. This is done to make several images (in 
     * time) of only one calculation, i.e. continue the calculation after saving 
     * the data.
     *
     * The integer MODE defines which scheme for the finite difference method 
     * should be used.
     */

    // calculate d_chi
    double d_chi = 1./((double)N_chi-1.);

    // some memory stuff, also create xvec.
    double * xvec = malloc(N_chi*sizeof(double));
    double ** final = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        final[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        xvec[i] = (double)i/((double)N_chi-1);
    }

    // solve the pde N_tau times, each step: use the result as the initial 
    // condition in the next step.
    for (int t=0; t<N_tau; t++) {
        switch(MODE) {
            case 0:
                // use the actual scheme
                coupled_variable_single( params, d_tau, d_chi, initial, final, N_chi );
                break;
            // use schemes from schemes.c
            case 1:
                cvs_nochain_fb( params, d_tau, d_chi, initial, final, N_chi );
                break;
            case 2:
                cvs_chain_centered( params, d_tau, d_chi, initial, final, N_chi );
                break;
            case 3:
                cvs_nochain_centered( params, d_tau, d_chi, initial, final, N_chi );
                break;
            case 4:
                cvs_chain_fb_second( params, d_tau, d_chi, initial, final, N_chi );
                break;
            case 5:
                cvs_nochain_fb_second( params, d_tau, d_chi, initial, final, N_chi );
                break;
            default:
                // use the actual scheme
                coupled_variable_single( params, d_tau, d_chi, initial, final, N_chi );
                break;
        }
        for (int i=0; i<N_chi; i++) {
            for (int j=0; j<4; j++) {
                initial[j][i] = final[j][i];
            }
        }
    }

    // save data to file.
    save_array_to_file(initial, xvec, fname, N_chi);

    // free the vector final.
    for (int j=0; j<4; j++) {
        free(final[j]);
    }
    free(final);
}

void coupled_stoichio( double * params, int N_chi, double d_tau, int N_tau,
        double ** initial, char * fname, double (*rfun)(double, double) ) {
    /*
     * params the same as in coupled_stoichio_single(), only an initial value 
     * for the configuration is needed. This is done to make several images (in 
     * time) of only one calculation, i.e. continue the calculation after saving 
     * the data.
     * As well as in coupled_stoichio_single(), the function pointer to the 
     * normalized function R(a,b) is needed `*rfun`.
     */

    // calculate d_chi
    double d_chi = 1./((double)N_chi-1.);

    // some memory stuff, also create xvec.
    double * xvec = malloc(N_chi*sizeof(double));
    double ** final = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        final[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        xvec[i] = (double)i/((double)N_chi-1);
    }

    // solve the pde N_tau times, each step: use the result as the initial 
    // condition in the next step.
    for (int t=0; t<N_tau; t++) {
        coupled_stoichio_single( params, d_tau, d_chi, initial, final,
                N_chi, rfun );
        for (int i=0; i<N_chi; i++) {
            for (int j=0; j<4; j++) {
                initial[j][i] = final[j][i];
            }
        }
    }

    // save data to file.
    save_array_to_file(initial, xvec, fname, N_chi);

    // free the vector final.
    for (int j=0; j<4; j++) {
        free(final[j]);
    }
    free(final);
}

void run_tests( char * prefix, double * Dfac, double * bdy, int N_tau,
        int tausteps ) {
    /*
     * a wrapper written to run and save the required tests if the solution of 
     * the pde is correct.
     * prefix is what the filenames get appended to
     * Dfac is an array to set D_b = Dfac[0]*D_a and D_c = Dfac[1]*D_a.
     * bdy is an array to set a_0 = bdy[0], a_1 = bdy[1]
     * N_tau is the number of maximal time steps taken
     * tausteps is the number of intermediary times that are saved.
     */

    // define the parameters. Mostly as in task d), but with some freedoms.
    double L = 1.;
    double D_a = 4.e-7;
    double D_b = Dfac[0]*D_a;
    double D_c = Dfac[1]*D_a;
    double R = 1.;
    double N_1 = R/10.;
    double N_2 = R/10.;
    double a_0 = bdy[0];
    double b_0 = bdy[1];
    double c_0 = 3.e-2;
    double t_0 = 1.;
    int N_chi = 2001;
    double d_tau = 0.1;

    // initialize the initial condition
    double ** initial = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        initial[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        for (int j=0; j<4; j++) {
            initial[j][i] = 0.0;
        }
    }
    // set params into array
    double params[8] = {t_0*D_a/pow(L,2), t_0*D_b/pow(L,2), t_0*D_c/pow(L,2),
        t_0*R*a_0, t_0*N_1*a_0, t_0*N_2*a_0, c_0/a_0, b_0/a_0};

    // run the thing.
    int currtau = 0;
    int taudelta = N_tau / tausteps;
    for (int T=0; T<tausteps; T++) {
        char fname[64];
        currtau += taudelta;
        sprintf(fname,"%sDfac0-%.3f_Dfac1-%.3f_a0-%.3f_b0-%.3f_t-%.3e.dat",
                prefix,Dfac[0],Dfac[1],a_0,b_0,(double)currtau*t_0);
        coupled_constant( params, N_chi, d_tau, taudelta, initial, fname );
    }

    // free the initial condition
    for (int j=0; j<4; j++) {
        free(initial[j]);
    }
    free(initial);
}

void run_constant ( char * prefix, int N_tau, int tausteps, double c_0,
        int N_chi ) {
    /*
     * a wrapper written to perform some tasks involving the simulation for 
     * constant D
     * prefix is what the filenames get appended to
     * N_tau is the number of maximal time steps taken
     * tausteps is the number of intermediary times that are saved.
     * c_0 is a parameter
     */

    // define the parameters. As in task d), only c_0 is free (and N_chi for 
    // checks), for later task e).
    double L = 1.; double D_a = 4.e-7; double D_b = 2./3.*D_a;
    double D_c = 8./15.*D_a; double R = 1.; double N_1 = R/10.;
    double N_2 = R/10.; double a_0 = 1.; double b_0 = 10.*a_0;
    double t_0 = 1.;

    double d_tau = 0.1;

    // initialize the initial condition
    double ** initial = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        initial[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        for (int j=0; j<4; j++) {
            initial[j][i] = 0.0;
        }
    }
    // set params into array
    double params[8] = {t_0*D_a/pow(L,2), t_0*D_b/pow(L,2), t_0*D_c/pow(L,2),
        t_0*R*a_0, t_0*N_1*a_0, t_0*N_2*a_0, c_0/a_0, b_0/a_0};

    // run the thing.
    int currtau = 0;
    int taudelta = N_tau / tausteps;
    for (int T=0; T<tausteps; T++) {
        char fname[64];
        currtau += taudelta;
        sprintf(fname,"%sc0-%.3e_t-%.3e.dat",
                prefix,c_0,(double)currtau*t_0);
        coupled_constant( params, N_chi, d_tau, taudelta, initial, fname );
    }

    // free the initial condition
    for (int j=0; j<4; j++) {
        free(initial[j]);
    }
    free(initial);
}

void run_variable ( char * prefix, int N_tau, int tausteps, double s_0,
        int N_chi, int MODE ) {
    /*
     * a wrapper written to perform some tasks involving the simulation for 
     * variable D
     * prefix is what the filenames get appended to
     * N_tau is the number of maximal time steps taken
     * tausteps is the number of intermediary times that are saved.
     * s_0 is a parameter
     *
     * The integer ´MODE ´defines which scheme to use; case statement for the 
     * choice of schemes in the function ´coupled_variable()´.
     */

    // define the parameters. As in task d), only c_0 is free (and N_chi for 
    // checks), for later task e).
    double L = 1.; double D_a = 4.e-7; double D_b = 2./3.*D_a;
    double D_c = 8./15.*D_a; double R = 1.; double N_1 = R/10.;
    double N_2 = R/10.; double a_0 = 1.; double b_0 = 10.*a_0;
    double t_0 = 1.; double c_0 = 3.e-2;

    double d_tau = 0.1;

    // initialize the initial condition
    double ** initial = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        initial[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        for (int j=0; j<4; j++) {
            initial[j][i] = 0.0;
        }
    }
    // set params into array
    double params[9] = {t_0*D_a/pow(L,2), t_0*D_b/pow(L,2), t_0*D_c/pow(L,2),
        t_0*R*a_0, t_0*N_1*a_0, t_0*N_2*a_0, c_0/a_0, b_0/a_0, s_0/a_0};

    // run the thing.
    int currtau = 0;
    int taudelta = N_tau / tausteps;
    for (int T=0; T<tausteps; T++) {
        char fname[64];
        currtau += taudelta;
        sprintf(fname,"%ss0-%.3e_t-%.3e.dat",
                prefix,s_0,(double)currtau*t_0);
        coupled_variable( params, N_chi, d_tau, taudelta, initial, fname, MODE );
    }

    // free the initial condition
    for (int j=0; j<4; j++) {
        free(initial[j]);
    }
    free(initial);
}

void run_stoichio( char * prefix, int N_tau, int tausteps, int N_chi,
        double st_a, double st_b, double (*rfun)(double, double),
        char * rfun_name ) {
    /*
     * a wrapper written to perform some tasks involving the simulation for 
     * different stoichiometries
     * prefix is what the filenames get appended to
     * N_tau is the number of maximal time steps taken
     * tausteps is the number of intermediary times that are saved.
     * rfun is a function pointer used for the reaction rate function, st_a and 
     * st_b are the respective stoichiometry constants: st_a a + st_b b -> c.
     */

    // define the parameters. As in task d), only c_0 is free (and N_chi for 
    // checks), for later task e).
    double L = 1.; double D_a = 4.e-7; double D_b = 2./3.*D_a;
    double D_c = 8./15.*D_a; double R = 1.; double N_1 = R/10.;
    double N_2 = R/10.; double a_0 = 1.; double b_0 = 10.*a_0;
    double t_0 = 1.; double c_0 = 3.e-2;

    double d_tau = 0.1;

    // initialize the initial condition
    double ** initial = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        initial[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        for (int j=0; j<4; j++) {
            initial[j][i] = 0.0;
        }
    }
    // set params into array
    double params[10] = {t_0*D_a/pow(L,2), t_0*D_b/pow(L,2), t_0*D_c/pow(L,2),
        t_0*R*a_0, t_0*N_1*a_0, t_0*N_2*a_0, c_0/a_0, b_0/a_0, st_a, st_b};

    // run the thing.
    int currtau = 0;
    int taudelta = N_tau / tausteps;
    for (int T=0; T<tausteps; T++) {
        char fname[64];
        currtau += taudelta;
        sprintf(fname,"%s%s_st-a-%i_st-b-%i_t-%.3e.dat",
                prefix, rfun_name, (int)st_a, (int)st_b ,(double)currtau*t_0);
        coupled_stoichio( params, N_chi, d_tau, taudelta, initial, fname, rfun );
    }

    // free the initial condition
    for (int j=0; j<4; j++) {
        free(initial[j]);
    }
    free(initial);
}

void run_constant_different_D ( char * prefix, int N_tau, int tausteps,
        int N_chi, double * params ) {
    /*
     * a wrapper written to perform simulations for task f)
     * prefix is what the filenames get appended to
     * N_tau is the number of maximal time steps taken
     * tausteps is the number of intermediary times that are saved.
     * params = {D_a_fac,D_b_fac,D_c_fac,a_0,b_0}
     */

    // define the parameters..
    double L = 1.; double R = 1.; double N_1 = R/10.; double N_2 = R/10.;
    double t_0 = 1.; double c_0 = 3.e-2;

    double D_a = 4.e-7*params[0];
    double D_b = 4.e-7*params[1];
    double D_c = 4.e-7*params[2];
    double a_0 = params[3]; double b_0 = params[4];

    double d_tau = 0.1;

    // initialize the initial condition
    double ** initial = malloc(4*sizeof(double*));
    for (int j=0; j<4; j++) {
        initial[j] = malloc(N_chi*sizeof(double));
    }
    for (int i=0; i<N_chi; i++) {
        for (int j=0; j<4; j++) {
            initial[j][i] = 0.0;
        }
    }
    // set params into array
    double sim_params[8] = {t_0*D_a/pow(L,2), t_0*D_b/pow(L,2), t_0*D_c/pow(L,2),
        t_0*R*a_0, t_0*N_1*a_0, t_0*N_2*a_0, c_0/a_0, b_0/a_0};

    // run the thing.
    int currtau = 0;
    int taudelta = N_tau / tausteps;
    for (int T=0; T<tausteps; T++) {
        char fname[64];
        currtau += taudelta;
        sprintf( fname, "%sDa-%.1f_Db-%.1f_Dc-%.1f_a0-%.1f_b0-%.1f_t-%.3e.dat",
                prefix, params[0], params[1], params[2], params[3], params[4],
                (double)currtau*t_0 );
        coupled_constant( sim_params, N_chi, d_tau, taudelta, initial, fname );
    }

    // free the initial condition
    for (int j=0; j<4; j++) {
        free(initial[j]);
    }
    free(initial);
}

int main( int argc, char ** argv ) {

    if (false) {
        // change to true for generating the data; just a check if the peaks are 
        // the same for higher N_chi.  They are.
        int N[5] = {201,401,801,1601,2001};
#pragma omp parallel for
        for (int i=0; i<5; i++) {
            char prefix[32];
            sprintf(prefix,"N_sim/N-%i",N[i]);
            run_constant(prefix, 2000000, 1, 3.e-2, N[i]);
        }
    }

    if (false) {
        // change to true if the simulation for constant D (different c0) should 
        // be run.
        double c0[3] = {3.e-2,1.e-2,6.e-2};
#pragma omp parallel for
        for (int i=0; i<3; i++) {
            run_constant( "const_sim/", 2000000, 20, c0[i], 2001 );
        }
    }

    if (false) {
        // change to true if the tests should be run.
        double Dfac[2] = {2.0,1.0};
        double bdy[2] = {1.0,1.0};
        run_tests( "tests/", Dfac, bdy, 100000, 10 );
        Dfac[0] = 0.5; bdy[1] = 0.5;
        run_tests( "tests/", Dfac, bdy, 100000, 10 );
    }

    if (false) {
        // change to true if the simulation for variable D should be run for a 
        // variety of different schemes all defined in schemes.c, and for 
        // different s0. Only used for testing purposes.
#define s0num 5
#define schemenum 6
        double s0[s0num] = {1.,3.5,5.1,7.1,10.1};
        int modes[schemenum] = {0, 1, 2, 3, 4, 5};
        char names[schemenum][128] = {
            "var_sim_schemes/chain-fb_",
            "var_sim_schemes/nochain-fb_",
            "var_sim_schemes/chain-cen_",
            "var_sim_schemes/nochain-cen_",
            "var_sim_schemes/chain-fb-sec_",
            "var_sim_schemes/nochain-fb-sec_"
        };
#pragma omp parallel for
        for (int i=0; i<schemenum*s0num; i++) {
            int k=i%schemenum; int l=i/schemenum;
            run_variable( names[k], 2000000, 5, s0[l], 501 , modes[k]);
        }
    }

    if (false) {
        // change to true if the simulation for variable D should be run (for 
        // different s0) for *just one scheme*, the one that was finally 
        // applied.
        double s0[5] = {100.,10.,1.,0.1,0.01};
#pragma omp parallel for
        for (int i=0; i<5; i++) {
            run_variable( "var_sim/", 2000000, 5, s0[i], 501, 0);
        }
    }

    if (false) {
        // change to true if the simulation for variable stoichiometries should 
        // be run. To make sure every case gets treated, set the environment 
        // variable OMP_NUM_THREADS to at least 6 before running.
        double R_one (double a, double b) {
            return a*b;
        }
        double R_two (double a, double b) {
            return pow(a,2)*b;
        }
#pragma omp parallel
        { int num = omp_get_thread_num();
            if (num == 0)
                run_stoichio( "st_sim/", 2000000, 5, 501, 2., 1.,
                        &R_one, "simple" );
            else if (num == 1)
                run_stoichio( "st_sim/", 2000000, 5, 501, 10., 1.,
                        &R_one, "simple" );
            else if (num == 2)
                run_stoichio( "st_sim/", 2000000, 5, 501, 1., 10.,
                        &R_one, "simple" );
            else if (num == 3)
                run_stoichio( "st_sim/", 2000000, 5, 501, 2., 1.,
                        &R_two, "powers" );
            else if (num == 4)
                run_stoichio( "st_sim/", 2000000, 5, 501, 10., 1.,
                        &R_two, "powers" );
            else if (num == 5)
                run_stoichio( "st_sim/", 2000000, 5, 501, 1., 10.,
                        &R_two, "powers" );
        }
    }

    if (false) {
        // change to true if some values for a set of different values of $D$ 
        // and $a_0$, $b_0$ should be generated.
        double D_b_vec[5] = {0.3,1.,3.,10.,20.};
        double D_c_vec[6] = {0.0,0.3,1.,3.,10.,20.};
        double b_0_vec[5] = {0.3,1.,3.,10.,20.};
#pragma omp parallel for
        for (int i=0; i<5; i++) {
            for (int j=0; j<6; j++) {
                for (int k=0; k<5; k++) {
                    double parameters[5] = {1.0, D_b_vec[i], D_c_vec[j], 3.0,
                        b_0_vec[k]};
                    run_constant_different_D( "D_sim/", 2000000, 1, 501,
                            parameters );
                }
            }
        }
    }
}
