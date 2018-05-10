/*
 * This file contains several schemes for testing the pde with variable D. They 
 * are not commented since they basically all do the same thing as the function 
 * ´coupled_variable_single()´ in the file reaction.c and only differ in the 
 * mathematical expressions for the numerical scheme itself
 */

// some prototypes needed.
double D_func( double, double );
double D_func_d( double, double );
double heaviside( double );

void cvs_nochain_fb ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    // forward-backward-averaged scheme without applying the chainrule
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    double sigma_0 = params[8];

    for (int i=0; i<N; i++) {
        int ip = shift_idx(i,1,N);
        int im = shift_idx(i,-1,N);
        double D_g = D_func( in[3][i], sigma_0 );
        double D_gp = D_func( in[3][ip], sigma_0 );
        double D_gm = D_func( in[3][im], sigma_0 );

        out[0][i] = in[0][i] + D_g * d_a*( in[0][shift_idx(i,1,N)] +
                in[0][shift_idx(i,-1,N)] - 2.*in[0][i] ) - r*in[0][i]*in[1][i] +
                0.5*d_a*( (D_g-D_gm) * (in[0][shift_idx(i,1,N)] - in[0][i]) + 
                          (D_gp-D_g) * (in[0][i] - in[0][shift_idx(i,-1,N)]) );
        out[1][i] = in[1][i] + D_g * d_b*( in[1][shift_idx(i,1,N)] +
                in[1][shift_idx(i,-1,N)] - 2.*in[1][i] ) - r*in[0][i]*in[1][i] +
                0.5*d_b*( (D_g-D_gm) * (in[1][shift_idx(i,1,N)] - in[1][i]) + 
                          (D_gp-D_g) * (in[1][i] - in[1][shift_idx(i,-1,N)]) );
        out[2][i] = in[2][i] + D_g * d_c*( in[2][shift_idx(i,1,N)] +
                in[2][shift_idx(i,-1,N)] - 2.*in[2][i] ) + r*in[0][i]*in[1][i] -
                n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i] +
                0.5*d_c*( (D_g-D_gm) * (in[2][shift_idx(i,1,N)] - in[2][i]) + 
                          (D_gp-D_g) * (in[2][i] - in[2][shift_idx(i,-1,N)]) );
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

void cvs_chain_centered ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    // centered scheme with applying the chainrule
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    double sigma_0 = params[8];

    for (int i=0; i<N; i++) {
        int ip = shift_idx(i,1,N);
        int im = shift_idx(i,-1,N);
        double D_g = D_func( in[3][i], sigma_0 );
        double D_g_d = D_func_d( in[3][i], sigma_0 );

        out[0][i] = in[0][i] + D_g * d_a*(in[0][ip] + in[0][im] - 2.*in[0][i]) -
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*D_g_d*(in[3][ip]-in[3][im])*(in[0][ip]-in[0][im]);
        out[1][i] = in[1][i] + D_g * d_b*(in[1][ip] + in[1][im] - 2.*in[1][i]) -
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*D_g_d*(in[3][ip]-in[3][im])*(in[1][ip]-in[1][im]);
        out[2][i] = in[2][i] + D_g * d_c*(in[2][ip] + in[2][im] - 2.*in[2][i]) +
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*D_g_d*(in[3][ip]-in[3][im])*(in[2][ip]-in[2][im])
                - n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

void cvs_nochain_centered ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    // centered scheme without applying the chainrule
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    double sigma_0 = params[8];

    for (int i=0; i<N; i++) {
        int ip = shift_idx(i,1,N);
        int im = shift_idx(i,-1,N);
        double D_g = D_func( in[3][i], sigma_0 );
        double D_gp = D_func( in[3][ip], sigma_0 );
        double D_gm = D_func( in[3][im], sigma_0 );

        out[0][i] = in[0][i] + D_g * d_a*(in[0][ip] + in[0][im] - 2.*in[0][i]) -
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*(D_gp-D_gm)*(in[0][ip]-in[0][im]);
        out[1][i] = in[1][i] + D_g * d_b*(in[1][ip] + in[1][im] - 2.*in[1][i]) -
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*(D_gp-D_gm)*(in[1][ip]-in[1][im]);
        out[2][i] = in[2][i] + D_g * d_c*(in[2][ip] + in[2][im] - 2.*in[2][i]) +
                    r*in[0][i]*in[1][i] +
                    0.25*d_a*(D_gp-D_gm)*(in[2][ip]-in[2][im])
                - n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

int shift_two_idx(int idx, int shift, int N) {
    int new_idx = idx + shift;
    if (new_idx < 0) {
        return 0;
    } else if (new_idx >= N) {
        return N;
    } else {
        return new_idx;
    }
}

void cvs_nochain_fb_second ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    // forward-backward-averaged scheme of second order in d_chi without 
    // applying the chainrule
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    double sigma_0 = params[8];

    for (int i=0; i<N; i++) {
        int ip = shift_two_idx(i,1,N);
        int im = shift_two_idx(i,-1,N);
        int ipp = shift_two_idx(i,2,N);
        int imm = shift_two_idx(i,-2,N);
        double D_g = D_func( in[3][i], sigma_0 );
        double D_gp = D_func( in[3][ip], sigma_0 );
        double D_gm = D_func( in[3][im], sigma_0 );
        double D_gpp = D_func( in[3][ipp], sigma_0 );
        double D_gmm = D_func( in[3][imm], sigma_0 );

        //double D_g_d = D_func_d( in[3][i], sigma_0 )/6./6./2.;

        out[0][i] = in[0][i] + D_g * d_a/12.*(
    - in[0][imm] + 16.*in[0][im] - 30.*in[0][i] + 16.*in[0][ip] - in[0][ipp] ) -
                    r*in[0][i]*in[1][i] +
                    d_a/6./6./2.*(
                        (D_gmm-6.*D_gm+3.*D_g+2.*D_gp) *
                        (-2.*in[0][im]-3.*in[0][i]+6.*in[0][ip]-in[0][ipp]) +
                        (-2.*D_gm-3.*D_g+6.*D_gp-D_gpp) *
                        (in[0][imm]-6.*in[0][im]+3.*in[0][i]+2.*in[0][ip])
                            );
        out[1][i] = in[1][i] + D_g * d_b/12.*(
    - in[1][imm] + 16.*in[1][im] - 30.*in[1][i] + 16.*in[1][ip] - in[1][ipp] ) -
                    r*in[0][i]*in[1][i] +
                    d_b/6./6./2.*(
                        (D_gmm-6.*D_gm+3.*D_g+2.*D_gp) *
                        (-2.*in[1][im]-3.*in[1][i]+6.*in[1][ip]-in[1][ipp]) +
                        (-2.*D_gm-3.*D_g+6.*D_gp-D_gpp) *
                        (in[1][imm]-6.*in[1][im]+3.*in[1][i]+2.*in[1][ip])
                            );
        out[2][i] = in[2][i] + D_g * d_c/12.*(
    - in[2][imm] + 16.*in[2][im] - 30.*in[2][i] + 16.*in[2][ip] - in[2][ipp] ) +
                    r*in[0][i]*in[1][i] +
                    d_c/6./6./2.*(
                        (D_gmm-6.*D_gm+3.*D_g+2.*D_gp) *
                        (-2.*in[2][im]-3.*in[2][i]+6.*in[2][ip]-in[2][ipp]) +
                        (-2.*D_gm-3.*D_g+6.*D_gp-D_gpp) *
                        (in[2][imm]-6.*in[2][im]+3.*in[2][i]+2.*in[2][ip])
                            )
                - n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}

void cvs_chain_fb_second ( double * params, double d_tau, double d_chi,
        double ** in, double ** out, int N ) {
    // forward-backward-averaged scheme of second order in d_chi with applying 
    // the chainrule
    double d_a = params[0] * d_tau / pow(d_chi,2);
    double d_b = params[1] * d_tau / pow(d_chi,2);
    double d_c = params[2] * d_tau / pow(d_chi,2);
    double r = params[3]*d_tau;
    double n_1 = params[4]*d_tau;
    double n_2 = params[5]*d_tau;

    double sigma_0 = params[8];

    for (int i=0; i<N; i++) {
        int ip = shift_two_idx(i,1,N);
        int im = shift_two_idx(i,-1,N);
        int ipp = shift_two_idx(i,2,N);
        int imm = shift_two_idx(i,-2,N);
        double D_g = D_func( in[3][i], sigma_0 );

        double D_g_d = D_func_d( in[3][i], sigma_0 )/6./6./2.;

        out[0][i] = in[0][i] + D_g * d_a/12.*(
    - in[0][imm] + 16.*in[0][im] - 30.*in[0][i] + 16.*in[0][ip] - in[0][ipp] ) -
                    r*in[0][i]*in[1][i] +
                    d_a*D_g_d*(
                        (in[3][imm]-6.*in[3][im]+3.*in[3][i]+2.*in[3][ip]) *
                        (-2.*in[0][im]-3.*in[0][i]+6.*in[0][ip]-in[0][ipp]) +
                        (-2.*in[3][im]-3.*in[3][i]+6.*in[3][ip]-in[3][ipp]) *
                        (in[0][imm]-6.*in[0][im]+3.*in[0][i]+2.*in[0][ip])
                            );
        out[1][i] = in[1][i] + D_g * d_b/12.*(
    - in[1][imm] + 16.*in[1][im] - 30.*in[1][i] + 16.*in[1][ip] - in[1][ipp] ) -
                    r*in[0][i]*in[1][i] +
                    d_b*D_g_d*(
                        (in[3][imm]-6.*in[3][im]+3.*in[3][i]+2.*in[3][ip]) *
                        (-2.*in[1][im]-3.*in[1][i]+6.*in[1][ip]-in[1][ipp]) +
                        (-2.*in[3][im]-3.*in[3][i]+6.*in[3][ip]-in[3][ipp]) *
                        (in[1][imm]-6.*in[1][im]+3.*in[1][i]+2.*in[1][ip])
                            );
        out[2][i] = in[2][i] + D_g * d_c/12.*(
    - in[2][imm] + 16.*in[2][im] - 30.*in[2][i] + 16.*in[2][ip] - in[2][ipp] ) +
                    r*in[0][i]*in[1][i] +
                    d_c*D_g_d*(
                        (in[3][imm]-6.*in[3][im]+3.*in[3][i]+2.*in[3][ip]) *
                        (-2.*in[2][im]-3.*in[2][i]+6.*in[2][ip]-in[2][ipp]) +
                        (-2.*in[3][im]-3.*in[3][i]+6.*in[3][ip]-in[3][ipp]) *
                        (in[2][imm]-6.*in[2][im]+3.*in[2][i]+2.*in[2][ip])
                            )
                - n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2) -
                n_2*in[2][i]*in[3][i];
        out[3][i] = in[3][i] + n_1*heaviside(in[2][i]-params[6])*pow(in[2][i],2)
                + n_2*in[2][i]*in[3][i];
    }
    out[0][0] = 1.; out[0][N-1] = 0.;
    out[1][0] = 0.; out[1][N-1] = params[7];
    out[2][0] = 0.; out[2][N-1] = 0.;
    out[3][0] = 0.; out[3][N-1] = 0.;
}
