Monte Carlo simulated annealing (MCSA) method was used to obtain the nonlinear fits to the ion-solvent (EC or PC) dimer binding energies from quantum chemistry calculations. 

The parameters are initially set to zero and then randomly sampled. For a better fit, the root-mean square error (rmse) between the quantum mechanics and classical binding energies was minimized by variation of the parameters. The selection of the parameters is based on the Metropolis algorithm. As the initial temperature is decreased (annealing), random changes in the energy are accepted with a Boltzmann probability condition. 

mcsa_fitting.py finds the fitted parameters. Running the script from terminal will look like this:

python mcsa_fitting.py --ion "ion name" --charge "ion charge" --j "min ion displacement" --nstep 'number of Monte Carlo steps' --tempr 'initial temperature'

plot.py plots the graphs of quantum and molecular mechanics (before&after fitting) energies.