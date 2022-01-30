import math

# Electrostatic interaction parameter set up
# E (kcal/mol) =  332 * q1*q2/r where q1 and q2 are charges and r is the distance between the molecules. 
ESC = 627.51 # 1 Hartree 627.51 kcal/mol
ANGAU = 1 / 0.0529177 # 1 bohr is 0.529177249 Angstrom
KB = 0.001987 # Boltzmann constant in kcal/mol

qc3 = 0.0954
qos = -0.4049
qc = 0.848701
qo = -0.548501
qh1 = 0.0797

# Intermolecular interaction (LJ, VdW energies) parameter set up
lmbda = 12.74

# sigma's are in Angstrom and epsilon's are in kcal/mol
sig_c3 = 0.34
sig_os = 0.3
sig_c = 0.34
sig_o = 0.296
sig_h1 = 0.2471

eps_c3 = 0.4577
eps_os = 0.7113
eps_c = 0.3598
eps_o = 0.87864
eps_h1 = 0.0657
