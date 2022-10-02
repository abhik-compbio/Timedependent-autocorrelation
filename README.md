# About
* This code can be used to calculate time dependent autocorrelation function for protein dihedral angles.
* The code is capable to calculate TDCF for both main chain and side chain dihedral angle of protein.

# Licensing

If you use the code for published works, please cite as
* Spatio-temporal coordination among functional residues in protein, Dutta et al., Scientific reports, 2017.
* Correlated dipolar and dihedral fluctuations in a protein, Moulick et al., Chemical Physics Letters, 2022.

# Required input file
We need 2 input file two run the code.
* "dyncorr_respair.dat" - this file contain the name of two output file and residue number of which we need to calculate timedependent auto correlation functions.
* "input-dataset.dat"- dihedral angle value of each dihedral angle obtained from molecular dynamics simulations.

# Compilation
To run the code:
* gfortan auto-res_dih_dyncorr.f90 -o a.out
* ./a.out

