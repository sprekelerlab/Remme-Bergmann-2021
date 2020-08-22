#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

double*** Make3DDoubleArray(int arraySizeX, int arraySizeY, int arraySizeZ) {
    double*** theArray;
    int i, j;
    theArray = (double***) mxCalloc(arraySizeX,sizeof(double**));
    if (theArray == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");    
    for (i = 0; i < arraySizeX; i++) {
        theArray[i] = (double**) mxCalloc(arraySizeY,sizeof(double*));
        if (theArray[i] == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");
        for (j = 0; j < arraySizeY; j++) {
            theArray[i][j] = (double*) mxCalloc(arraySizeZ,sizeof(double));
            if (theArray[i][j] == NULL) mexErrMsgTxt("Can not allocate temporary variables\n");
        }
    }
    return theArray;
}

gsl_rng * Rng;               
const gsl_rng_type * Rng_T;

double uniform() {
    return gsl_rng_uniform(Rng);
}

unsigned long int uniform_int(unsigned long int n) {
	return gsl_rng_uniform_int(Rng,n);
}

double gaussian(double sigma) {
	return gsl_ran_gaussian(Rng,sigma);
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    double  *W_HPC_init, *W_init, *params;
    int     nsec, N_layer, N_cell, it_del, seed;
    double  dt, A_mean, A_std, w_max, lambda,lambda_q, alpha, tau_LTP, tau_LTD;

	if (nrhs != 3) mexErrMsgTxt("Error in number of input arguments!\n");
	if (nlhs != 2) mexErrMsgTxt("Error in number of output variables!\n");

    W_HPC_init    	= mxGetPr(prhs[0]);
    W_init    		= mxGetPr(prhs[1]);
    params    		= mxGetPr(prhs[2]);	
    nsec            = (int) params[0]; 
    dt              = params[1];
	N_layer 		= (int) params[2];
	N_cell			= (int) params[3];
	A_mean			= params[4];
	A_std			= params[5];
	it_del 			= (int) params[6];
	w_max 			= params[7];
	lambda			= params[8];
	lambda_q 		= params[9];
	alpha 			= params[10];
	tau_LTP			= params[11];
	tau_LTD 		= params[12];
	seed 			= (long) params[13];
    
    gsl_rng_env_setup();
    Rng_T 	= gsl_rng_default;
    Rng 	= gsl_rng_alloc(Rng_T);   
    gsl_rng_set(Rng, seed);

    double** 	A_in 	= Make2DDoubleArray(N_cell,N_layer);
    double** 	A_out 	= Make2DDoubleArray(N_cell,N_layer);
    double*  	A_HPC 	= Make1DDoubleArray(N_cell);	    		 	/* activity of HPC cells through W_HPC matrix */        
    double** 	W_HPC 	= Make2DDoubleArray(N_cell,N_cell);    		 	/* weight matrix from input layer 0 to HPC */        
    double*** 	W 		= Make3DDoubleArray(N_cell,N_cell,N_layer);           

    double*** 	A_in_hist  = Make3DDoubleArray(N_cell,N_layer,it_del); 	/* history of all input layer cells */
    double*** 	A_out_hist = Make3DDoubleArray(N_cell,N_layer,it_del); 	/* history of all output layer cells */
    double**  	A_HPC_hist = Make2DDoubleArray(N_cell,it_del);   		/* history of HPC cells */
    
    double* 	d_LTP	= Make1DDoubleArray(N_layer);
    double* 	d_LTD	= Make1DDoubleArray(N_layer);
	
    int  		Nind	= 1/dt;

    mwSignedIndex dims[3];
    dims[0] = N_cell;
    dims[1] = N_cell;
    dims[2] = N_layer; 
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N_cell, N_layer+2, mxREAL);

    double  	*W_save, *A_save;
    W_save   	= mxGetPr(plhs[0]);
	A_save 		= mxGetPr(plhs[1]);
		
    /* INITIALIZE */
    int i = 0;
    int j = 0;
	int k = 0;
	int hist_store_ind = 0;
	int hist_use_ind = 0;

	/* initialize learning rates */
	for (k=0;k<N_layer;k++) { 				/* Note numbering: Layer 1 is the fastest */
		d_LTP[k] = lambda*pow(lambda_q,k);
		/*d_LTP[k] = lambda;*/
		d_LTD[k] = alpha*d_LTP[k];
	}

    /* initialize weight matrices */
    for (i=0;i<N_cell;i++) {
        for (j=0;j<N_cell;j++) {
            W_HPC[i][j] = W_HPC_init[j+N_cell*i];
        }
    }
    for (k=0;k<N_layer;k++) {
	    for (i=0;i<N_cell;i++) {
	        for (j=0;j<N_cell;j++) {
	            W[i][j][k] = W_init[j+N_cell*i+k*N_cell*N_cell];
	        }
	    }
	}
	    
  	double** LTP_A = Make2DDoubleArray(N_cell,N_layer);  /* low pass filtered activities of input layer cells */
    double** LTD_A = Make2DDoubleArray(N_cell,N_layer);  /* low pass filtered activities of output layer cells */
    
    int 	sec  	= 0;
    int 	ind     = 0;
 		
	/* SECONDS LOOP */    
    for (sec=1;sec<=nsec;sec++) {
		/*mexPrintf("%d\n",sec);*/

        /* DT LOOP */
        for (ind=0;ind<Nind;ind++) {
			hist_store_ind = ind%it_del;
			hist_use_ind   = (ind+1)%it_del;	

            /* activity of 'sensory' input layer */
			k = N_layer - 1;
            for (i=0;i<N_cell;i++) {
                A_in[i][k] = A_mean + gaussian(A_std);
				if (A_in[i][k]<0) A_in[i][k] = 0;
			 	A_in_hist[i][k][hist_store_ind] = A_in[i][k]; /* store A_in activity */
            }
			
			/* activities of subsequent input layers (from periphery to HPC) */
			for (k=N_layer-2;k>=0;k--) {
	            for (i=0;i<N_cell;i++) {
	                A_in[i][k] = A_in_hist[i][k+1][hist_use_ind]; /* set activity using previous act in previous layer */
				 	A_in_hist[i][k][hist_store_ind] = A_in[i][k]; /* store A_in activity */
	            }			
			}

			/* activity of HPC cells */
            for (i=0;i<N_cell;i++) {
                A_HPC[i] = 0;
                for (j=0;j<N_cell;j++) {
                    A_HPC[i] += W_HPC[i][j]*A_in_hist[j][0][hist_use_ind];  /* from input layer closest to HPC, 1 delay ago */
                }
			 	A_HPC_hist[i][hist_store_ind] = A_HPC[i]; /* store A_HPC activity */
            }
					
            /* activities of first output layer (closest to HPC) */
            k = 0;
            for (i=0;i<N_cell;i++) {
                A_out[i][k] = 0;
                for (j=0;j<N_cell;j++) {
                    A_out[i][k] += 0.5*W[i][j][k]*A_in_hist[j][k][hist_use_ind]; /* the k=0 is input layer 1 (the closest to HPC) */
                }                
                A_out[i][k] 	+= 0.5*A_HPC_hist[i][hist_use_ind]; /* identity input of delayed HPC activity */
			 	A_out_hist[i][k][hist_store_ind] = A_out[i][k]; /* store A_out activity */
           	}
			
			/* activities of subsequent output layers (from HPC to periphery) */
			for (k=1;k<N_layer;k++) { 	/* layer 0 is already set above */
	            for (i=0;i<N_cell;i++) {
	                A_out[i][k] = 0;
	                for (j=0;j<N_cell;j++) {
	                    A_out[i][k]	+= 0.5*W[i][j][k]*A_in_hist[j][k][hist_use_ind]; /* note that input layers are counted from the HPC */
	                }                
	                A_out[i][k] 	+= 0.5*A_out_hist[i][k-1][hist_use_ind]; /* identity using delayed previous output layer activities */
				 	A_out_hist[i][k][hist_store_ind] = A_out[i][k]; /* store A_out activity */
	            }			
			}

            /* STDP of shortcut fibers */
			for (k=0;k<N_layer;k++) {
	            for (i=0;i<N_cell;i++) {
	                LTP_A[i][k] += dt*(A_in_hist[i][k][hist_use_ind] - LTP_A[i][k])/tau_LTP; /* this should represent low pass filtered activity at the synapse terminal, hence uses the propagated activity*/
	                LTD_A[i][k] += dt*(A_out[i][k] - LTD_A[i][k])/tau_LTD; /* this is low pass filtered activity at the post-synaptic site, and hence uses current activity */
	            }
			}
			for (k=0;k<N_layer;k++) {
	            for (i=0;i<N_cell;i++) {
	                for (j=0;j<N_cell;j++) {
	                    W[i][j][k] += dt*( A_out[i][k]*LTP_A[j][k]*d_LTP[k] - LTD_A[i][k]*A_in_hist[j][k][hist_use_ind]*d_LTD[k] );
	                    if (W[i][j][k]>w_max) W[i][j][k] = w_max;                    
	                    if (W[i][j][k]<0) W[i][j][k] = 0;                                  
	                }
	            }
			}		
			
        } /* end of DT loop */
    } /* end of sec loop */
	
	/* save shortcut weight matrices */
	for (k=0;k<N_layer;k++) {
	    for (i=0;i<N_cell;i++) {
	        for (j=0;j<N_cell;j++) {
	            W_save[j+i*N_cell+k*N_cell*N_cell]  = W[i][j][k];
	        }
	    }
	}
	/* save layer activities from final iteration */
    for (i=0;i<N_cell;i++) {
        A_save[i]= A_in[i][0];
    }
    for (i=0;i<N_cell;i++) {
        A_save[i+N_cell] = A_HPC[i];
    }	
    for (k=0;k<N_layer;k++) {
        for (i=0;i<N_cell;i++) {
            A_save[i+(k+2)*N_cell] = A_out[i][k];
        }
    }
}
