# Randomization of Approximate Bilinear Computation for Matrix Multiplication
This repository provides the code we use in our preprint paper 

O. A. Malik and S. Becker. Randomization of Approximate Bilinear Computation for Matrix Multiplication. *arXiv:1905.07439 [cs.DS]*, 17 May 2019.

It is available at https://arxiv.org/abs/1905.07439

## Paper abstract
We present a method for randomizing a formula for bilinear computation of matrix products. We consider the implications of such randomization when the formula itself is approximate, and when the formula is exact but its computation is plagued by numerical error due to finite precision arithmetic. Our theoretical results and numerical experiments indicate that our method can improve performance in both settings for a negligible increase in computational complexity.

## Requirements
Parts of our code relies on the following software:
* Tensor Toolbox version 2.6 by Bader, Kolda and others (available at http://www.sandia.gov/~tgkolda/TensorToolbox/).
* Round with significant digits by François Beauducel (available at https://www.mathworks.com/matlabcentral/fileexchange/26212-round-with-significant-digits).

## Installation
In addition to the required software listed above, this repository contains some C code that needs to be compiled. This can be done by executing **compile_rand_mat_mult_C.m** from the Matlab command line.

## Experiments
There are a total of six script files in this repository that we used to run the various experiments in our paper. We now describe which script was used for which experiment.
* **experiment1a.m** and **experiment1b.m**: Used to do compute some of the quantities mentioned in Section 2.2.
* **experiment2.m**: Used for the first experiment in Section 3.1 (Figure 1).
* **experiment3.m**: Used for the second experiment in Section 3.1 (Figure 2) and the first set of experiments in Section B of the supplement (Figures 6-9).
* **experiment4.m**: Used for the first experiment in Section 3.2 (Figure 3) and the second set of experiments in Section B of the supplement (Figures 10-13).
* **experiment5.m**: Used for the second experiment in Section 3.2 (Figure 4 and 5) and the third set of experiments in Section B of the supplement (Figures 14-17).

## Referencing this code
If you use our code in any of your own work, please provide a reference to this repository and/or our paper, as appropriate.

We have done our best to include relevant references in our code to make it easier to find original resources.

## Author contact information
Please feel free to contact me at any time if you have any questions or would like to provide feedback on this code or on our paper. I can be reached at osman.malik@colorado.edu.

## Licenses
This code uses the function roundsd by François Beauducel, which is available on MathWorks File Exchange. The license of that software is available in its original form in the folder random-approximate-matrix-multiplication/roundsd.

All other code in this project falls under the license in the root of this project.