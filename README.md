# Timedependent-autocorrelation
* This code can be used to calculate time dependent autocorrelation function for protein dihedral angles.
* The code is capable to calculate TDCF for both main chain and side chain dihedral angle of protein.
We need 2 input file two run the code.
* "dyncorr_respair.dat" - this file contain the name of two output file and residue number of which we need to calculate timedependent auto correlation functions.
* "input-dataset.dat"- dihedral angle value of each dihedral angle obtained from molecular dynamics simulations.

To run the code:

gfortan auto-res_dih_dyncorr.f90 -o a.out
./a.out

