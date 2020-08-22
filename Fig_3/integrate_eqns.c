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

double* Make1DDoubleArray(int arraySizeX) {
    double* theArray;
    int i;
    theArray = (double*) mxCalloc(arraySizeX,sizeof(double));
    if (theArray == NULL) mexErrMsgTxt("Cannot allocate temporary variables\n");
    return theArray;
}

double** Make2DDoubleArray(int arraySizeX, int arraySizeY) {
    double** theArray;
    int i;
    theArray = (double**) mxCalloc(arraySizeX,sizeof(double*));
    if (theArray == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");    
    for (i = 0; i < arraySizeX; i++) {
        theArray[i] = (double*) mxCalloc(arraySizeY,sizeof(double));
        if (theArray[i] == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");    
    }
    return theArray;
}

gsl_rng *Rng;               
const gsl_rng_type * Rng_T;

double uniform() {
    return gsl_rng_uniform(Rng);
}

unsigned long int uniform_int(unsigned long int n) {
	return gsl_rng_uniform_int(Rng,n);
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

    double  *W_SC_po_init, *W_PP_po_init, *W_PP_pp_init, *W_PPS_po_init, *W_PPS_pp_init, *params, *CA3_p_pf, *CA3_o_A, *EC_p_gf, *EC_o_A;
    int     nsec, N, seed, N_obj, it_del_CA3, it_del_CA1;
    double  dt, W_PP_max, W_PPS_max, W_SC_id, W_SCS_id, lambda_PP,lambda_PPS, alpha, tau_LTP, tau_LTD, env_len, CA3_pf_var, A_EC_max, A_CA3_max;

	if (nrhs != 10) mexErrMsgTxt("Error in number of input arguments!\n");
	if (nlhs != 4) mexErrMsgTxt("Error in number of output variables!\n");

    W_SC_po_init    = mxGetPr(prhs[0]);
    W_PP_po_init    = mxGetPr(prhs[1]);
    W_PP_pp_init    = mxGetPr(prhs[2]);	
    W_PPS_po_init   = mxGetPr(prhs[3]);
    W_PPS_pp_init   = mxGetPr(prhs[4]);	
    params          = mxGetPr(prhs[5]);
    CA3_p_pf        = mxGetPr(prhs[6]);
    CA3_o_A         = mxGetPr(prhs[7]);
    EC_p_gf         = mxGetPr(prhs[8]);
    EC_o_A          = mxGetPr(prhs[9]);
    nsec            = (int) params[0]; 
    dt              = params[1]; 
    N               = (int) params[2];
	W_PP_max 		= params[3];
	W_PPS_max 		= params[4];
    W_SC_id         = params[5];
	W_SCS_id 		= params[6];
    it_del_CA3      = (int) params[7];
	it_del_CA1 		= (int) params[8];
	lambda_PP       = params[9];
	lambda_PPS      = params[10];
    alpha           = params[11];
    tau_LTP         = params[12];
    tau_LTD         = params[13];
    env_len         = params[14];
    CA3_pf_var      = params[15];
    A_EC_max        = params[16];
    A_CA3_max       = params[17];
	N_obj 			= params[18];	
    seed            = (long) params[19];
    
    gsl_rng_env_setup();
    Rng_T = gsl_rng_default;
    Rng = gsl_rng_alloc(Rng_T);    
    gsl_rng_set (Rng, seed);

    double*  A_EC_p  = Make1DDoubleArray(N);                
    double*  A_EC_o  = Make1DDoubleArray(N);                 
    double*  A_CA3_p = Make1DDoubleArray(N);                
    double*  A_CA3_o = Make1DDoubleArray(N);             
    double*  A_CA1_p = Make1DDoubleArray(N);             
    double*  A_SUB_p = Make1DDoubleArray(N);             
    
    double** W_PP_pp = Make2DDoubleArray(N, N);           
    double** W_PP_po = Make2DDoubleArray(N, N);            
    double** W_SC_po = Make2DDoubleArray(N, N);            
    double** W_PPS_pp = Make2DDoubleArray(N, N);           
    double** W_PPS_po = Make2DDoubleArray(N, N);            

    double** del_CA3_p = Make2DDoubleArray(it_del_CA3, 2);     
    int*  	 del_CA3_o = Make1DIntArray(it_del_CA3);     
    double** del_CA1_p = Make2DDoubleArray(it_del_CA1, N);     
    
    double  d_LTP_PP    = lambda_PP;
    double  d_LTD_PP    = alpha*lambda_PP;
    double  d_LTP_PPS   = lambda_PPS;
    double  d_LTD_PPS   = alpha*lambda_PPS;
    
    int     Nind        = 1/dt;

    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, N, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, N, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(N, N, mxREAL);

    double  *W_PP_po_save, *W_PPS_po_save, *W_PP_pp_save, *W_PPS_pp_save;
    W_PP_po_save   	= mxGetPr(plhs[0]);
    W_PPS_po_save  	= mxGetPr(plhs[1]);
    W_PP_pp_save  	= mxGetPr(plhs[2]);
    W_PPS_pp_save  	= mxGetPr(plhs[3]);
    
    /* INITIALIZE */
    int i = 0;
    int j = 0;

    /* initialize weight matrices */
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            W_SC_po[j][i] = W_SC_po_init[j+N*i];
            W_PP_pp[j][i] = W_PP_pp_init[j+N*i];
            W_PP_po[j][i] = W_PP_po_init[j+N*i];
            W_PPS_pp[j][i] = W_PPS_pp_init[j+N*i];
            W_PPS_po[j][i] = W_PPS_po_init[j+N*i];
        }
    }
    
    double* LTD_PP_p 	= Make1DDoubleArray(N);
  	double* LTP_PP_p  	= Make1DDoubleArray(N);
    double* LTP_PP_o  	= Make1DDoubleArray(N);
    double* LTD_PPS_p  	= Make1DDoubleArray(N);
    
    double LTP_ud   = exp(-dt/tau_LTP);
    double LTD_ud   = exp(-dt/tau_LTD);
    
    int sec     = 0;
    int ind     = 0;
    
    double phi  = 0;
    double r    = 0;
    double xp   = 0;
    double yp   = 0;
    double xp_del = 0;
    double yp_del = 0;
    int o_ind	  = 0;
    int o_ind_del = 0;

        
    /* store initial W_PP_po + W_PPS_po matrices */
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            W_PP_po_save[j+i*N+0]  = W_PP_po[j][i];
            W_PPS_po_save[j+i*N+0] = W_PPS_po[j][i];
        }
    }
    
	/* SECONDS LOOP */    
    for (sec=1;sec<=nsec;sec++) {
        
        /* DT LOOP */
        for (ind=0;ind<Nind;ind++) {
                        
            /* generate position on circle */
            phi = 2*M_PI*uniform();
            r   = uniform() + uniform();    // TODO: why not the same formula as in matlab code? (sqrt(uniform)...)
            if (r>1) r = 2-r;
            xp  = env_len/2*r*cos(phi);
            yp  = env_len/2*r*sin(phi);
            del_CA3_p[ind%it_del_CA3][0] = xp;
            del_CA3_p[ind%it_del_CA3][1] = yp;           
            xp_del = del_CA3_p[(ind+1)%it_del_CA3][0];
            yp_del = del_CA3_p[(ind+1)%it_del_CA3][1];
            
            /* generate object */
			o_ind = uniform_int(N_obj);
            del_CA3_o[ind%it_del_CA3] = o_ind;
            o_ind_del = del_CA3_o[(ind+1)%it_del_CA3];
            
            /* activity of EC cells */
            for (i=0;i<N;i++) {
                A_EC_p[i] = 0;
                A_EC_p[i] =  cos((cos(           EC_p_gf[i+2*N])*(xp-EC_p_gf[i])  
                                + sin(           EC_p_gf[i+2*N])*(yp-EC_p_gf[i+N]))*EC_p_gf[i+3*N]);
                A_EC_p[i] += cos((cos(M_PI/3+    EC_p_gf[i+2*N])*(xp-EC_p_gf[i])  
                                + sin(M_PI/3+    EC_p_gf[i+2*N])*(yp-EC_p_gf[i+N]))*EC_p_gf[i+3*N]);
                A_EC_p[i] += cos((cos(2*M_PI/3+  EC_p_gf[i+2*N])*(xp-EC_p_gf[i])  
                                + sin(2*M_PI/3+  EC_p_gf[i+2*N])*(yp-EC_p_gf[i+N]))*EC_p_gf[i+3*N]);
                A_EC_p[i] = A_EC_max*(A_EC_p[i]+1.5)/4.5;
                        
				A_EC_o[i] = EC_o_A[i+N*o_ind];
            }

            /* activity of CA3 cells */
            for (i=0;i<N;i++) {
                A_CA3_p[i] = A_CA3_max*exp(-((xp_del-CA3_p_pf[i])*(xp_del-CA3_p_pf[i])+(yp_del-CA3_p_pf[i+N])*(yp_del-CA3_p_pf[i+N]))/(2*CA3_pf_var));
				A_CA3_o[i] = CA3_o_A[i+N*o_ind_del];				
            }

            /* activity of CA1 cells */
            for (i=0;i<N;i++) {
                A_CA1_p[i] = 0;
                for (j=0;j<N;j++) {
                    A_CA1_p[i] += W_PP_pp[i][j]*A_EC_p[j];
                    A_CA1_p[i] += W_PP_po[i][j]*A_EC_o[j];
                    A_CA1_p[i] += W_SC_po[i][j]*A_CA3_o[j];
                }                
                A_CA1_p[i] += W_SC_id*A_CA3_p[i]; /* identity */                
            }
			
			/* store CA1_p activity to implement delay to SUB */
			for (i=0;i<N;i++) del_CA1_p[ind%it_del_CA1][i] = A_CA1_p[i];
			
            /* activity of SUB cells */
            for (i=0;i<N;i++) {
                A_SUB_p[i] = 0;
                for (j=0;j<N;j++) {
                    A_SUB_p[i] += W_PPS_pp[i][j]*A_EC_p[j];
                    A_SUB_p[i] += W_PPS_po[i][j]*A_EC_o[j];
                }                
                A_SUB_p[i] += W_SCS_id*del_CA1_p[(ind+1)%it_del_CA1][i]; /* identity using delayed CA1 activity */                
            }			
            
            /* STDP of PP fibers */
            for (i=0;i<N;i++) {
                LTP_PP_p[i] += dt*(A_EC_p[i] - LTP_PP_p[i])/tau_LTP;
                LTP_PP_o[i] += dt*(A_EC_o[i] - LTP_PP_o[i])/tau_LTP;
            }
            for (i=0;i<N;i++) {
                LTD_PP_p[i] += dt*(A_CA1_p[i] - LTD_PP_p[i])/tau_LTD;
                LTD_PPS_p[i]+= dt*(A_SUB_p[i] - LTD_PPS_p[i])/tau_LTD;
            }
									
            for (i=0;i<N;i++) {
                for (j=0;j<N;j++) {
                    W_PP_pp[j][i] += dt*(A_CA1_p[j]*LTP_PP_p[i]*d_LTP_PP - LTD_PP_p[j]*A_EC_p[i]*d_LTD_PP);
                    W_PP_po[j][i] += dt*(A_CA1_p[j]*LTP_PP_o[i]*d_LTP_PP - LTD_PP_p[j]*A_EC_o[i]*d_LTD_PP);
                    if (W_PP_pp[j][i]>W_PP_max) W_PP_pp[j][i] = W_PP_max;                    
                    if (W_PP_po[j][i]>W_PP_max) W_PP_po[j][i] = W_PP_max;
                    if (W_PP_pp[j][i]<0) W_PP_pp[j][i] = 0;                                  
                    if (W_PP_po[j][i]<0) W_PP_po[j][i] = 0;
                    W_PPS_pp[j][i] += dt*(A_SUB_p[j]*LTP_PP_p[i]*d_LTP_PPS - LTD_PPS_p[j]*A_EC_p[i]*d_LTD_PPS);
                    W_PPS_po[j][i] += dt*(A_SUB_p[j]*LTP_PP_o[i]*d_LTP_PPS - LTD_PPS_p[j]*A_EC_o[i]*d_LTD_PPS);
                    if (W_PPS_pp[j][i]>W_PPS_max) W_PPS_pp[j][i] = W_PPS_max;                    
                    if (W_PPS_po[j][i]>W_PPS_max) W_PPS_po[j][i] = W_PPS_max;
                    if (W_PPS_pp[j][i]<0) W_PPS_pp[j][i] = 0;                                  
                    if (W_PPS_po[j][i]<0) W_PPS_po[j][i] = 0;
                }
            }
        }
    }
	/* store weight matrices */
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            W_PP_po_save[j+i*N]  = W_PP_po[j][i];
            W_PPS_po_save[j+i*N] = W_PPS_po[j][i];			
            W_PP_pp_save[j+i*N]  = W_PP_pp[j][i];
            W_PPS_pp_save[j+i*N] = W_PPS_pp[j][i];
        }
    } 	
}
