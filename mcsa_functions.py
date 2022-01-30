# functions used by the mcsa fitting script

import math
from mcsa_parameters import *


def vdw_pot(i, a_j, b_j, c_j, p, t):
    x, y, z, n = t
    a_i, b_i, c_i = p
    # xn, yn and zn are the coordinates of the ion
    xn = x[n - 1]
    yn = y[n - 1]
    zn = z[n - 1]

    xr = xn - x[i]
    yr = yn - y[i]
    zr = zn - z[i]

    b = 2 / (1 / math.pow(b_j, 6) + 1 / math.pow(b_i, 6))
    b_ij = math.pow(b, 0.16666)
    a_ij = math.sqrt(a_j * a_i) * b / math.pow(b_j, 3) / math.pow(b_i, 3)
    c_ij = math.sqrt(c_j * c_i)
    d = 0.00005 * 4.184
    aa = d * math.pow(12.0 / b_ij, 12)

    rnL = math.sqrt(xr * xr + yr * yr + zr * zr)

    return a_ij * math.exp(-1.0 * b_ij * rnL) - (c_ij / math.pow(rnL, 6)) + aa / math.pow(rnL, 12)


def a_par(eps, lmbda):
    return 6 * eps * math.exp(lmbda) / (lmbda - 6)


def b_par(sig, lmbda):
    return lmbda / sig


def c_par(sig, eps, lmbda):
    return eps * lmbda * math.pow(sig, 6) / (lmbda - 6)


def calc_coulomb(qi, t):

    def el_c(i, qi, q_atom, t):
        x, y, z, n = t
        # xn, yn and zn are the coordinates of the ion
        xn = x[n - 1]
        yn = y[n - 1]
        zn = z[n - 1]

        xr = (xn - x[i]) * ANGAU
        yr = (yn - y[i]) * ANGAU
        zr = (zn - z[i]) * ANGAU
        rnL = math.sqrt(xr * xr + yr * yr + zr * zr)

        return qi * q_atom / rnL
        
    elc = 0
    elc = el_c(0, qi, qc3, t)
    elc = elc + el_c(1, qi, qc3, t)
    elc = elc + el_c(2, qi, qos, t)
    elc = elc + el_c(3, qi, qc, t)
    elc = elc + el_c(4, qi, qo, t)
    elc = elc + el_c(5, qi, qos, t)
    elc = elc + el_c(6, qi, qh1, t)
    elc = elc + el_c(7, qi, qh1, t)
    elc = elc + el_c(8, qi, qh1, t)
    elc = elc + el_c(9, qi, qh1, t)

    elc = ESC * elc

    return elc


def calc_vdw(p_i, t):
    # Buckingham potential parameters
    a_c3 = a_par(eps_c3, lmbda)
    a_os = a_par(eps_os, lmbda)
    a_c = a_par(eps_c, lmbda)
    a_o = a_par(eps_o, lmbda)
    a_h1 = a_par(eps_h1, lmbda)

    b_c3 = b_par(sig_c3, lmbda)
    b_os = b_par(sig_os, lmbda)
    b_c = b_par(sig_c, lmbda)
    b_o = b_par(sig_o, lmbda)
    b_h1 = b_par(sig_h1, lmbda)

    c_c3 = c_par(sig_c3, eps_c3, lmbda)
    c_os = c_par(sig_os, eps_os, lmbda)
    c_c = c_par(sig_c, eps_c, lmbda)
    c_o = c_par(sig_o, eps_o, lmbda)
    c_h1 = c_par(sig_h1, eps_h1, lmbda)

    vdw = 0
    vdw = vdw + vdw_pot(0, a_c3, b_c3, c_c3, p_i, t)
    vdw = vdw + vdw_pot(1, a_c3, b_c3, c_c3, p_i, t)
    vdw = vdw + vdw_pot(2, a_os, b_os, c_os, p_i, t)
    vdw = vdw + vdw_pot(3, a_c, b_c, c_c, p_i, t)
    vdw = vdw + vdw_pot(4, a_o, b_o, c_o, p_i, t)
    vdw = vdw + vdw_pot(5, a_os, b_os, c_os, p_i, t)
    vdw = vdw + vdw_pot(6, a_h1, b_h1, c_h1, p_i, t)
    vdw = vdw + vdw_pot(7, a_h1, b_h1, c_h1, p_i, t)
    vdw = vdw + vdw_pot(8, a_h1, b_h1, c_h1, p_i, t)
    vdw = vdw + vdw_pot(9, a_h1, b_h1, c_h1, p_i, t)

    # Multiplying by 0.23901 converts vdW potential to kcal/mol
    return 0.23901 * vdw
