from math import pi
from dataclasses import dataclass
from _utils import viscosity_bloodplasma


@dataclass
class FlowModelParameters:
    # r_c_bar = 0.001/2.       # Quarteroni, Table 1.7 bore/2
    r_c_bar = 3.75e-4          #
    l_c_bar = 0.02             # Quarteroni, Table 1.7

    n_rev = 460

    rho_c = 997 * 1e-3
    rho_t = 1060 * 1e-3

    n_sev = n_rev / 3 # every capillary leaves at two ends somewhere

    K_c = n_sev * r_c_bar**4 * pi / 8 / (1e-2)
    # K_c = 1e-13 * 1e4  # open?
    K_t = 1e-18 * 1e4  # Tobias paper homogenisierung


    mu_c = viscosity_bloodplasma(r_c_bar)       # y
    mu_t = 0.79722 * 1e-3 * 10                  # ?

    S_ct = 460 * r_c_bar * l_c_bar * 2 * pi / 1e-3  # fromula from Josza

    #L_cv = 1e-8  # guess
    #L_cv = 5e-9  # guess
    #L_cv = 2.5e-9  # guess
    #L_cv = 1e-8  # guess
    L_cv = 5e-9  # guess
    L_ct = 1e-9  # guess
    #L_tl = 1e-8  # guess
    L_tl = 5e-9  # guess

    sigma = 0.  # ?
    pi_int = 6.6e3
    pi_bl = 3.3e4

    p_ven = 10 * 1333
    p_lym = 1333.2

@dataclass
class FlowModelParametersOld:
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
    D_c = 1.7e-5                # siehe unten -
    D_t = 1.7e-5                # d'Angleo -> D_{t,b} cm^2/s

    lambda_t = 6 * 1e-8 * 1e3   # d'Angelo theta