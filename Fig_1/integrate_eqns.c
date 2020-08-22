#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

int* Make1DIntArray(int arraySizeX) {
    int* theArray;
    int i;
    theArray = (int*) mxCalloc(arraySizeX,sizeof(int));
    if (theArray == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");
    return theArray;
}

int** Make2DIntArray(int arraySizeX, int arraySizeY) {
    int** theArray;
    int i;
    theArray = (int**) mxCalloc(arraySizeX,sizeof(int*));
    if (theArray == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");

    for (i = 0; i < arraySizeX; i++) {
        theArray[i] = (int*) mxCalloc(arraySizeY,sizeof(int));
        if (theArray[i] == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");
    }
    return theArray;
}

double* Make1DDoubleArray(int arraySizeX) {
    double* theArray;
    int i;
    theArray = (double*) mxCalloc(arraySizeX,sizeof(double));
    if (theArray == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");
    return theArray;
}

gsl_rng *Rng;
const gsl_rng_type * Rng_T;
int     gsl_ran_poisson(const gsl_rng * r,double mu);
void    gsl_ran_choose (const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size);
float   getpid();

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    double  *we_sc_init, *we_pp_init, *params;
    double  *we_sc_result, *we_pp_result, *rate_result, *v_result;
    double  v, ge_sc, ge_pp;
    int     nsec, Ne, save_we_ivl;
    double  dt, E_freq, we_max, delay_sc_pp, lambda_sc, lambda_pp, alpha;


    if (nrhs != 3) mexErrMsgTxt("Error in number of input arguments!\n");
    if (nlhs != 4) mexErrMsgTxt("Error in number of output variables!\n");

    we_sc_init      = mxGetPr(prhs[0]);
    we_pp_init      = mxGetPr(prhs[1]);
    params          = mxGetPr(prhs[2]);
    nsec            = (int) params[0];
    dt              = params[1];
    Ne              = params[2];
    E_freq          = params[3];
    we_max          = params[4];
    delay_sc_pp     = params[5];
    lambda_sc       = params[6];
    lambda_pp       = params[7];
    alpha           = params[8];
    save_we_ivl     = (int) params[9];

    gsl_rng_env_setup();
    Rng_T = gsl_rng_default;
    Rng = gsl_rng_alloc(Rng_T);
    long seed;

    seed = time (NULL) * getpid();
    gsl_rng_set (Rng, seed);

    /* some hard-coded parameters */
    double  taum    = 20;
    double  El      = -70;
    double  Ee      = 0;
    double  Vt      = -54; /* threshold voltage */
    double  Vr      = -60; /* reset potential */
    double  tref    = 1.75; /* total duration of refractory period */
    double  tspike  = 0.8; /* duration of depolarized phase of spike; for plotting purposes */
    double  tau_ge  = 5; /* synaptic decay time constant */
    double  tau_LTP = 20;
    double  tau_LTD = 20;

    double  d_LTP_sc   = lambda_sc;
    double  d_LTD_sc   = alpha*lambda_sc;
    double  d_LTP_pp   = lambda_pp;
    double  d_LTD_pp   = alpha*lambda_pp;

    /* initialize output arrays */
    plhs[0] = mxCreateDoubleMatrix(Ne, ceil(1.0*nsec/save_we_ivl), mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Ne, ceil(1.0*nsec/save_we_ivl), mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nsec, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,(int)(1000/dt), mxREAL);

    we_sc_result    = mxGetPr(plhs[0]);
    we_pp_result    = mxGetPr(plhs[1]);
    rate_result     = mxGetPr(plhs[2]);
    v_result        = mxGetPr(plhs[3]);

    v       = -70;
    ge_sc   = 0;
    ge_pp   = 0;

    /* update factors for exponentially decaying variables */
    double ge_ud    = exp(-dt/tau_ge);
    double LTP_ud   = exp(-dt/tau_LTP);
    double LTD_ud   = exp(-dt/tau_LTD);

    /* initialize variables */
    double* we_sc   = Make1DDoubleArray(Ne);
    double* we_pp   = Make1DDoubleArray(Ne);
    int i = 0;
    for (i=0;i<Ne;i++) {we_sc[i] = we_sc_init[i]; we_pp[i] = we_pp_init[i];}

    double* LTP_sc  = Make1DDoubleArray(Ne);
    double* LTP_pp  = Make1DDoubleArray(Ne);
    double LTD_sc   = 0;
    double LTD_pp   = 0;

    int refrac      = 0;
    double firing   = 0;
    double t_after_onset = 0;

    double dv1, dv2, v_tmp;

    int sec     = 0;
    double t    = 0;
    int ind     = 0;
    int Nind    = 1000/dt;
    int ind2    = 0;
    int rnd_tmp = 0;
    int ind_tmp = 0;
    int rate    = 0;
    int input_indcs[Ne];
    for (i=0;i<Ne; i++) {input_indcs[i]=i;}

    double E_prob    = Ne * dt/1000 * E_freq;
    int    delay_ind = delay_sc_pp/dt;

    int  matrix_wd      = ceil(10*E_prob);
    int* Nfired_e_sc_vec = Make1DIntArray(1000/dt);
    int* Nfired_e_pp_vec = Make1DIntArray(1000/dt);
    int** fired_e_sc_mat = Make2DIntArray(1000/dt,matrix_wd);
    int** fired_e_pp_mat = Make2DIntArray(1000/dt,matrix_wd);

    /* seconds loop */
    for (sec=0;sec<nsec;sec++) {

        /* copy delay_sc_pp milliseconds from previous second to start of next second */
        for (ind=0;ind<delay_ind;ind++){
            ind_tmp = Nind - delay_ind + ind;
            Nfired_e_sc_vec[ind] = Nfired_e_pp_vec[ind_tmp];
            for (ind2=0;ind2<Nfired_e_pp_vec[ind_tmp];ind2++){
                fired_e_sc_mat[ind][ind2] = fired_e_pp_mat[ind_tmp][ind2] ;
            }
        }

        /* generate pp input spiketimes for the entire second */
        for (ind=0;ind<Nind;ind++) {
            rnd_tmp = gsl_ran_poisson(Rng,E_prob); /* how many pp spikes in timestep ind */
            if (rnd_tmp>matrix_wd) rnd_tmp = matrix_wd;
            Nfired_e_pp_vec[ind] = rnd_tmp;
            gsl_ran_choose(Rng, fired_e_pp_mat[ind], rnd_tmp, input_indcs, Ne, sizeof(int)); /* which pp inputs are active in timestep ind */
        }

        /* copy pp input spiketimes to sc input with a delay of delay_sc_pp */
        for (ind=0;ind<(Nind - delay_ind);ind++){
            Nfired_e_sc_vec[ind + delay_ind] = Nfired_e_pp_vec[ind];
            for (ind2=0;ind2<Nfired_e_pp_vec[ind];ind2++){
                fired_e_sc_mat[ind + delay_ind][ind2] = fired_e_pp_mat[ind][ind2] ;
            }
        }

        rate = 0;
        /* dt loop*/
        for (ind=0;ind<Nind;ind++) {
            /* which sc inputs are active in timestep ind */
            for (i=0;i<Nfired_e_sc_vec[ind];i++) {
                LTP_sc[fired_e_sc_mat[ind][i]] += d_LTP_sc;
                we_sc[fired_e_sc_mat[ind][i]]  += LTD_sc;
                if (we_sc[fired_e_sc_mat[ind][i]]<0) we_sc[fired_e_sc_mat[ind][i]] = 0;
            }
            /* which pp inputs are active in timestep ind */
            for (i=0;i<Nfired_e_pp_vec[ind];i++) {
                LTP_pp[fired_e_pp_mat[ind][i]] += d_LTP_pp;
                we_pp[fired_e_pp_mat[ind][i]]  += LTD_pp;
                if (we_pp[fired_e_pp_mat[ind][i]]<0) we_pp[fired_e_pp_mat[ind][i]] = 0;
            }

            if (refrac) {
                t_after_onset = (t-firing);
                if (t_after_onset<tspike) {v = 30;} /* add spike for plotting purposes */
                else {v = Vr;}
                if (t_after_onset>tref) refrac = 0;
            } else if (v>Vt) {
                rate        +=1;
                v           = 30;  /* add spike for plotting purposes */
                refrac      = 1;
                firing      = t;

                LTD_sc     -= d_LTD_sc;
                for (i=0;i<Ne;i++) {
                    we_sc[i] += LTP_sc[i];
                    if (we_sc[i]>we_max) we_sc[i] = we_max;
                }

                LTD_pp     -= d_LTD_pp;
                for (i=0;i<Ne;i++) {
                    we_pp[i] += LTP_pp[i];
                    if (we_pp[i]>we_max) we_pp[i] = we_max;
                }
            }

            for (i=0;i<Nfired_e_sc_vec[ind];i++)  ge_sc += we_sc[fired_e_sc_mat[ind][i]];
            for (i=0;i<Nfired_e_pp_vec[ind];i++)  ge_pp += we_pp[fired_e_pp_mat[ind][i]];

            if (!refrac){ /* Heun's method */
                dv1     = ( El-v       + (ge_sc+ge_pp)*(Ee-v)      )/taum;
                v_tmp   = v + dt*dv1;
                dv2     = ( El-v_tmp   + (ge_sc+ge_pp)*(Ee-v_tmp)  )/taum;
                v      += dt*0.5*(dv1+dv2);
            }

            for (i=0;i<Ne;i++){
                LTP_sc[i]   *= LTP_ud;
                LTP_pp[i]   *= LTP_ud;
            }

            ge_sc       *= ge_ud;
            ge_pp       *= ge_ud;
            LTD_sc      *= LTD_ud;
            LTD_pp      *= LTD_ud;

            if(sec==(10)){
                v_result[ind]   = v;
            }
            t += dt;
        }

        if ((sec%10)==0) {
            mexPrintf("time = %d \t rate = %d \n",sec,rate);
            mexEvalString("pause(.0001);");
        }
        rate_result[sec] = rate;
        if ((sec%save_we_ivl)==0) {
            for (i=0;i<Ne;i++){
                we_sc_result[(int)(sec/save_we_ivl)*Ne+i]  = we_sc[i];
                we_pp_result[(int)(sec/save_we_ivl)*Ne+i]  = we_pp[i];
            }
        }
    }

    mxFree(we_sc);
    mxFree(we_pp);
    mxFree(LTP_sc);
    mxFree(LTP_pp);
}

