#*****************************************************************************
# Monte Carlo Simulated Annealing (MCSA) fitting developed by Ayse Arslanargin
# University of Cincinnati, OH, Department of Physics & Chemistry, 2014-2015
#*****************************************************************************

import argparse

import math
import random
from mcsa_functions import *
from mcsa_parameters import *


def main(params):
    ion = params.ion
    qi = int(params.qi)
    j = int(params.j)
    nstep = int(params.nstep)
    tempr0 = int(params.tempr0)

    n = 11

    rmse_best = 999
    rmse_old = 999

    sig_old = 0.0
    eps_old = 0.0

    a_best = 0
    b_best = 0
    c_best = 0

    sig = 0
    eps = 0
    coord_dict = {}


    FILES_ROOT = f'./ions/{ion}'

    # Reading quantum binding energies
    qm_file = open(FILES_ROOT + '-ec-qm.dat', 'r')
    qm_lines = qm_file.readlines()
    qm_ener = [float(line.split(",")[1]) for line in qm_lines]

    # Reading coordinates of ion-solvent pairs to calculate vdw&coulomb energies and putting them in a dictionary for easy access
    for s in range(-j, 40):
        c_x = []
        c_y = []
        c_z = [] 

        coord_file = open(FILES_ROOT + '/' + str(s) + '/output.gro', 'r')
        coord_lines = coord_file.readlines()
        coord_lines = coord_lines[2:-1]

        for line in coord_lines:
            coord = line.split()[3:]
            c_x.append(float(coord[0]))
            c_y.append(float(coord[1]))
            c_z.append(float(coord[2]))

            coord_dict[s] = (c_x, c_y, c_z)

    # Starting fitting algorithm
    for step in range(nstep):
        en_diff = 0
        rmse = 0
        tempr = tempr0 * math.exp(- 4 * step / nstep)
        
        weight = []
        x = []
        y = []
        z = []

        # Generate new trial points
        sig = sig_old + 0.0001 * random.random()
        eps = eps_old + 0.0001 * random.random()

        for s in range(-j, 40):
            x = coord_dict[s][0]
            y = coord_dict[s][1]
            z = coord_dict[s][2]
            t = (x, y, z, n)

            a_i = a_par(eps, lmbda)
            b_i = b_par(sig, lmbda)
            c_i = c_par(sig, eps, lmbda)
            p_i = (a_i, b_i, c_i)

            coulomb = calc_coulomb(qi, t)
            vdW = calc_vdw(p_i, t)

            mm_ener = coulomb + vdW
            # Boltzmann weighting function applied to all points for better fitting (optimizing the fit near the potential minimum)
            power = -1.0 * qm_ener[s+j] / (KB * tempr)
            # Required to prevent overflow 
            if power < 709:               
                weight.append(math.exp(power))
            else:
                return print(f'Parameters for {ion} are A = {p_best[0] / 4.18:.1f} B = {p_best[1] / 10:.3f} C = {p_best[2] / 4.18 * 10**6:.1f} Sigma = {sig_best: .2f} Epsilon = {eps_best: .2f}')

            en_diff = math.pow((mm_ener - qm_ener[s+j]), 2)

            rmse += weight[s+j] * en_diff

        rmse = math.sqrt(rmse / sum(weight))
        d_rmse = rmse - rmse_old

        if step == 0:
            # Required to prevent overflow
            p = 0   
        else:
            d_rmse = rmse - rmse_old
            boltz = -1 * rmse / (KB * tempr)
            p = math.exp(boltz)
        p0 = math.exp(-1 * d_rmse / (KB * tempr))

        if d_rmse < 0.0:
            p = 1.0
            accept = True
        elif p > p0:
            # Cost function (system energy) increases, so the chage is accepted with probability.
            accept = True
        else:
            accept = False

        if accept:
            rmse_old = rmse
            if rmse <= rmse_best:
                # Cost function (system energy) decreases, so change in the parameters is accepted.
                rmse_best = rmse
                p_best = (a_i, b_i, c_i)
                sig_best, eps_best = sig, eps
        else:
            # If the accept is rejected return to the previous accept
            sig, eps  = sig_old, eps_old 

        sig_old, eps_old = sig, eps 

        if (step % 500) == 0:
            print(f'{step} {tempr:.1f} {accept:1d} {rmse:.3f} {rmse_best:.3f} {p_best[0]:.6f} {p_best[1]:.2f} {p_best[2]:.5f} {sig_best:.2f} {eps_best:.2f}')
        
    print(f'Parameters for {ion} are A = {p_best[0] / 4.18:.1f} B = {p_best[1] / 10:.3f} C = {p_best[2] / 4.18 * 10**6:.1f} Sigma = {sig_best: .2f} Epsilon = {eps_best: .2f}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fitting parameters')

    parser.add_argument('--ion', help='ion name')
    parser.add_argument('--qi', help='ion charge')
    parser.add_argument('--j', help='min value of the ion displacement')
    parser.add_argument('--nstep', help='number of Monte Carlo steps')
    parser.add_argument('--tempr0', help='initial temperature')

    args = parser.parse_args()

    main(args)
