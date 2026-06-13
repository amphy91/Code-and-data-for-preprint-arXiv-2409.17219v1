# Code-and-data-for-Fermionic-parton-theory-of-Rydberg-Z-quantum-spin-liquid-
This repository contains the codes and data files for the paper "Fermionic parton theory of Rydberg Z₂ quantum spin liquid"

Contents:
1. The file "8a_dsf.f90" contains the code for calculating the data for Fig.3. The corresponding data file is given by "dsf.dat".
2. The file "8a_ssf.f90" contains the code for calculating the data for Fig.2(b). The corresponding data file is given by "ssf.dat".
3.There are two mathematica notebooks "DSF.nb" and "SSF.nb" which produce the corresponding plot with appropriate colorscheme.

Requirements:
The codes were developed using Gfortran and mathematica
Required packages:
blas and lapack

Reproducing the results:
1. Run "8a_dsf.f90" and "8a_ssf.f90" to generate the data file using the library blas and lapack.
2. Use the mathematica notebooks for plotting.

Contact:
For question regarding the code, please contact atanu.maity@uni-wuerzburg.de/amphy91@gmail.com

