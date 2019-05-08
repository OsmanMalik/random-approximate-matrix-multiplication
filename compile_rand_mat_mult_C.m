%COMPILE_RAND_MAT_MULT_C Compiles C function rand_mat_mult_C.c
%
%   COMPILE_RAND_MAT_MULT_C properly compiles the function
%   rand_mat_mult_C.c with the blas flag. It compiles the source code in
%   two different versions, one using floats and one using doubles. This is
%   done by compiling with and without the -DSINGLE flag. The outputs in
%   each case are named appropriately so that the compiled files work
%   together with the function rand_mat_mult_C_wrapper. Note that blas is
%   required for the compilation to work properly.

mex rand_mat_mult_C.c -lmwblas -DSINGLE -output rand_mat_mult_C_single
mex rand_mat_mult_C.c -lmwblas -output rand_mat_mult_C_double