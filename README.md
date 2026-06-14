# Code-and-data-for-Fermionic-parton-theory-of-Rydberg-Z-quantum-spin-liquid-
This repository contains the codes and data files for the paper "Fermionic parton theory of Rydberg Z₂ quantum spin liquid"

Contents:
1. The file "Ansatze_Table-I.nb" contains the mathematica notebook for calculating the symmetry allowed mean-field parameters given in Table I using the PSG solutions given in the supplimental information.
2. The file "RubyWFZ.jl" contains the Julia/ITensor implementation of the DMRG calculation for the Rydberg Hamiltonian on the ruby lattice. The calculation is performed on a 24 × 2 ruby-lattice cluster, corresponding to 288 sites, at parameters Δ/Ω = 1.7 and V/Ω = 50, with interactions included up to third-nearest neighbors. This code generates the DMRG ground state and evaluates the local magnetization and equal-time spin-spin correlation function in the S^z channel.
The file "Sz_24_2_D1.7_V50.0.txt" contains the site-resolved expectation values ⟨S_i^z⟩ obtained from the DMRG ground state. The file "ZZ_24_2_D1.7_V50.0.txt" contains the full real-space correlation matrix ⟨S_i^z S_j^z⟩. These DMRG data are used to compute the static structure factor shown in Fig. 2(a) of the manuscript.
3. The file "dsf.f90" contains the code for calculating the data for Fig.3. The corresponding data file is given by "dsf8a.dat".
4. The file "ssf.f90" contains the code for calculating the data for Fig.2(b). The corresponding data file is given by "ssf8a.dat".
4.There are two mathematica notebooks "DSF.nb" and "SSF.nb" which produce the corresponding plot with appropriate colorscheme.

Requirements:
The mean-field codes were developed using gfortran and Mathematica, with BLAS and LAPACK required for the Fortran calculations. The DMRG code was developed in Julia using the ITensor library, together with the packages MKL, ITensors, ITensors.HDF5, DelimitedFiles, Printf, and LinearAlgebra.

Reproducing the results:
1. Use "Ansatze_Table-I.nb" for reproducing the symmetry allowed mean-field parameters given in Table I.
2. julia RubyWFZ.jl 24 2 17 50
where the arguments specify Nx = 24, Ny = 2, Δ/Ω = 1.7, and V/Ω = 50.
3. Run "dsf.f90" and "ssf.f90" to generate the data file using the library blas and lapack.
5. Use the mathematica notebooks for plotting.

Contact:
For question regarding the code, please contact atanu.maity@uni-wuerzburg.de/amphy91@gmail.com

