/*
// RAND_MAT_MULT_C.C	Computes A*B approximately
//	
//	C = RAND_MAT_MULT_C(A, B, Y, epsilon, S, P) returns an approximation
//	of the product A*B computed via a bilinear computation defined by
//	the factor matrices in Y, the adjustment factor epsilon, the 
//	realization of Rademacher random variables in S, and the realization of
//	permutation random variables in P. 
//
//	If compiled with (without) -DSINGLE flag then:
//		- A and B should be matrices of type single (double), 
//		- Y should be a cell containing single (double) matrices, 
//		- epsilon should also be a single (double) precision scalar, 
//
//	S and P should always be cells containing matrices of type int32.
//
// COMPILATION:
//	To compile for single mode, run "mex rand_mat_mult_C.c -lmwblas -DSINGLE -output filename"
//	To compile for double mode, run "mex rand_mat_mult_C.c -lmwblas -output filename"
*/

#include "mex.h"
#include "blas.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Define variable type number to be either float or double */
#ifdef SINGLE
typedef float number;
#else
typedef double number;
#endif

/* Declare global variables for inputs */
number *A;
number *B;
number **Y;
number *epsilon;
int32_T **S;
int32_T **P;
mwIndex no_rec;
mwIndex n;
mwIndex R;

/* Declare global variables for outputs */
number *C;
	
/* Function for computing output C */
void compute_output(number *C_sub, number *A_sub, number *B_sub, size_t N_sub, mwIndex rec_level) {
	mwIndex r, c, i, j, k, l, m;
	number *A_block;
	number *B_block;
	number *C_block;
	number *C_block_temp;
	size_t block_sz;
	size_t block_numel;
	size_t unit_stride;
	number Y1_val, Y2_val, Y3_val;
	char *chn = "N";
	number alpha, beta;

	if(rec_level > no_rec) { /* Naive matrix multiply */
		/* Use xgemm from blas */
		alpha = 1.0;
		beta = 0.0;
#ifdef SINGLE
		sgemm(chn, chn, &N_sub, &N_sub, &N_sub, &alpha, A_sub, &N_sub, B_sub, &N_sub, &beta, C_sub, &N_sub);
#else
		dgemm(chn, chn, &N_sub, &N_sub, &N_sub, &alpha, A_sub, &N_sub, B_sub, &N_sub, &beta, C_sub, &N_sub);
#endif
		/* Old code
		for(c = 0; c < N_sub; ++c) {
			for(r = 0; r < N_sub; ++r) {
				C_sub[r + N_sub*c] = 0;
				for(i = 0; i < N_sub; ++i) {
					C_sub[r + N_sub*c] += A_sub[r + N_sub*i]*B_sub[i + N_sub*c];
				}
			}
		}
		*/
	} else { /* Fast algorithm */
		block_sz = N_sub/n;
		block_numel = block_sz*block_sz;
		unit_stride = 1;
		for(j = 0; j < n; ++j) {
			for(i = 0; i < n; ++i) {
				C_block = calloc(block_numel, sizeof(number));
				for(r = 0; r < R; ++r) {
					Y3_val = Y[2][(P[rec_level-1][0+3*i]-1)*n + (P[rec_level-1][2+3*j]-1) + n*n*r];
					if(Y3_val != 0) {
						A_block = calloc(block_numel, sizeof(number));
						B_block = calloc(block_numel, sizeof(number));
						for(k = 0; k < n; ++k) {
							for(l = 0; l < n; ++l) {
								Y1_val = Y[0][(P[rec_level-1][0+3*k]-1) + (P[rec_level-1][1+3*l]-1)*n + n*n*r];
								if(Y1_val != 0) {
									for(c = 0; c < block_sz; ++c) {
										for(m = 0; m < block_sz; ++m) {
											A_block[m + c*block_sz] += Y1_val * S[rec_level-1][0+3*k] * S[rec_level-1][1+3*l] * A_sub[m + k*block_sz + N_sub*(c + l*block_sz)];
										}
									}
								}
								Y2_val = Y[1][(P[rec_level-1][1+3*k]-1) + (P[rec_level-1][2+3*l]-1)*n + n*n*r];
								if(Y2_val != 0) {
									for(c = 0; c < block_sz; ++c) {
										for(m = 0; m < block_sz; ++m) {
											B_block[m + c*block_sz] += Y2_val * S[rec_level-1][1+3*k] * S[rec_level-1][2+3*l] * B_sub[m + k*block_sz + N_sub*(c + l*block_sz)];
										}
									}
								}
							}
						}
						C_block_temp = malloc(block_numel*sizeof(number));
						compute_output(C_block_temp, A_block, B_block, block_sz, rec_level+1);
						
						/* Use blas code for adding Y3_val*C_block_temp to C_block */
#ifdef SINGLE
						saxpy(&block_numel, &Y3_val, C_block_temp, &unit_stride, C_block, &unit_stride);
#else
						daxpy(&block_numel, &Y3_val, C_block_temp, &unit_stride, C_block, &unit_stride);
#endif
						
						/* Old code
						for(c = 0; c < N_sub/n; ++c) {
							for(m = 0; m < N_sub/n; ++m) {
								C_block[m + c*N_sub/n] += Y3_val*C_block_temp[m + c*N_sub/n];
							}
						}
						*/
						
						free(C_block_temp);
						free(B_block);
						free(A_block);
					}
				}
				for(l = 0; l < block_sz; ++l) {
					for(k = 0; k < block_sz; ++k) {
						C_sub[k + i*block_sz + N_sub*(l + j*block_sz)] = 1.0/(1.0+*epsilon) * S[rec_level-1][0+3*i] * S[rec_level-1][2+3*j] * C_block[k + block_sz*l];
					}
				}
				free(C_block);
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Declare other variables */
	mwIndex i;
	mwIndex no_cols;
	mwIndex no_rows;
	number C_result;
	number *mw_array;
	
	/* Get input variables */
	A = (number *) mxGetData(prhs[0]);
	B = (number *) mxGetData(prhs[1]);
	Y = malloc(3*sizeof(number *));
	for(i = 0; i < 3; ++i) {
		Y[i] = (number *) mxGetData(mxGetCell(prhs[2], i));
	}
	n = sqrt(mxGetM(mxGetCell(prhs[2], 0)));
	R = mxGetN(mxGetCell(prhs[2], 0));
	epsilon = (number *) mxGetData(prhs[3]);
	no_rec = mxGetNumberOfElements(prhs[4]);
	S = malloc(no_rec*sizeof(int32_T *));
	P = malloc(no_rec*sizeof(int32_T *));
	for(i = 0; i < no_rec; ++i) {
		S[i] = (int32_T *) mxGetData(mxGetCell(prhs[4], i));
		P[i] = (int32_T *) mxGetData(mxGetCell(prhs[5], i));
	}

	/* Create the output matrix */
	no_rows = mxGetM(prhs[0]);
	no_cols = mxGetN(prhs[0]);
#ifdef SINGLE
	plhs[0] = mxCreateNumericMatrix(no_rows, no_cols, mxSINGLE_CLASS, mxREAL);
#else
	plhs[0] = mxCreateDoubleMatrix(no_rows, no_cols, mxREAL);
#endif
	C = (number *) mxGetData(plhs[0]);
	
	/* Compute output */
	compute_output(C, A, B, no_rows, 1);
	
	/* Free dynamically allocated memory */
	free(P);
	free(S);
	free(Y);
}