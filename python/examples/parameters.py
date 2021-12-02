from math import pi
from dataclasses import dataclass


@dataclass
class FlowModelParameters:
    rho_c = 997 * 1e-3
    rho_t = 1060 * 1e-3
    K_c = 1e-13 * 1e4  # open?
    # K_c = 1e-12 * 1e4  # wildest guess!!!
    K_t = 1e-18 * 1e4  # Tobias paper homogenisierung
    mu_c = 0.0805438
    mu_t = 0.79722 * 1e-3 * 10

    # L_cv = 1e-7  # guess
    # L_cv = 0.667 / 1333  # Ottesen, Olufsen, Larsen, p. 153
    # L_cv = 4.641e-7  # Liverpool!
    L_cv = 1e-8  # guess
    #L_cv = 2.2e-3 * 10-4  # p.142 carlo

    p_ven = 10 * 1333

    # L_ct = 1e-6  # ?
    # L_ct = 1e-11  # -> L_cap Tobias paper
    L_ct = 1e-9  # guess

    r_c_bar = 3.375e-4
    l_c_bar = 0.06
    S_ct = 460 * r_c_bar * l_c_bar * 2 * pi / 1e-3  # ?

    sigma = 0.  # ?
    pi_int = 6.6e3
    pi_bl = 3.3e4

    L_tl = 1e-8  # ?
    p_lym = 1333.2


@dataclass
class TransportModelParameters:
    D_c = 1.7e-5    # siehe unten -
    D_t = 1.7e-5    # d'Angleo -> D_{t,b} cm^2/s
    lambda_t = 5 * 1e-8 * 1e3   # d'Angelo theta