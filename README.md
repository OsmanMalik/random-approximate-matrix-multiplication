# Randomization of Approximate Bilinear Computation for Matrix Multiplication
This repository provides the code we use in our paper 

Osman Asif Malik & Stephen Becker (2020) *Randomization of approximate bilinear computation for matrix multiplication*, International Journal of Computer Mathematics: Computer Systems Theory, DOI: [10.1080/23799927.2020.1861104](https://doi.org/10.1080/23799927.2020.1861104)

## Paper abstract
We present a method for randomizing formulas for bilinear computation of matrix products which does not increase the leading order complexity of the computation. We consider the implications of such randomization when there are two sources of error. The first source is due to the computation formula itself only being approximately correct. Such formulas come up when numerically searching for faster matrix multiplication algorithms. The second source is due to using floating point arithmetic. This kind of error is especially important when computing on low precision hardware like GPUs. Our theoretical results and numerical experiments indicate that our method can improve performance when the two kinds of error are present individually, as well as when they are present at the same time.

## Requirements
Parts of our code relies on the following software:
* Tensor Toolbox version 2.6 by Bader, Kolda and others (available at http://www.sandia.gov/~tgkolda/TensorToolbox/).
* Round with significant digits by François Beauducel (available at https://www.mathworks.com/matlabcentral/fileexchange/26212-round-with-significant-digits).

## Installation
In addition to the required software listed above, this repository contains some C code that needs to be compiled. This can be done by executing **compile_rand_mat_mult_C.m** from the Matlab command line.

## Experiments
There are a total of seven script files in this repository that we used to run the various experiments in our paper. We now describe which script was used for which experiment.
* **experiment1a.m** and **experiment1b.m**: Used to compute some of the quantities mentioned in Section 2.2.
* **experiment2.m**: Used for the first experiment in Section 3.1 (Figure 1).
* **experiment3.m**: Used for the second experiment in Section 3.1 (Figures 2-4) and the first set of experiments in the appendix (Figures A1-A5).
* **experiment4.m**: Used for the first experiment in Section 3.2 (Figures 5-6) and the second set of experiments in the appendix (Figures A6-A10).
* **experiment5.m**: Used for the second experiment in Section 3.2 (Figures 7-13) and the third set of experiments in the appendix (Figures A11-A15).
* **experiment6.m**: Used for the experiment in Section 3.3 (Figure 14).

## Referencing this code
If you use our code in any of your own work, please reference our paper:
```
@article{malik-becker-ijcm-cst,
  author = { Osman Asif   Malik  and  Stephen   Becker },
  title = {Randomization of approximate bilinear computation for matrix multiplication},
  journal = {International Journal of Computer Mathematics: Computer Systems Theory},
  pages = {1-40},
  year  = {2020},
  publisher = {Taylor & Francis},
  doi = {10.1080/23799927.2020.1861104},
  URL = {https://doi.org/10.1080/23799927.2020.1861104},
  eprint = {https://doi.org/10.1080/23799927.2020.1861104}
}
```

We have done our best to include relevant references in our code to make it easier to find original resources.

## Author contact information
Please feel free to contact me at any time if you have any questions or would like to provide feedback on this code or on our paper. I can be reached at osman.malik@colorado.edu.

## Licenses
This code uses the function roundsd by François Beauducel, which is available on MathWorks File Exchange. The license of that software is available in its original form in the folder random-approximate-matrix-multiplication/roundsd.

All other code in this project falls under the license in the root of this project.
