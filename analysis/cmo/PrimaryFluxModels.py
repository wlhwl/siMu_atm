import numpy as np


def GST3(Z, E: np.array):
    # dN/dE = \Sum a*E^(-gamma-1) * exp(-E / Z / R)
    flux_single_population = lambda AGammaR_List: AGammaR_List[0] * np.power(E,-AGammaR_List[1]-1) * np.exp(-E/(Z*AGammaR_List[2]))
    pars = []
    if Z==1:
        pars = [
            # a, gamma, R
            [7e3, 1.661, 1.2e5], # Pop 1
            [150, 1.4, 4e6],     # Pop 2
            [1.4, 1.4, 1.3e9],   # Pop 3
        ]
    elif Z==2:
        pars = [
            # a, gamma, R
            [3200, 1.58, 1.2e5], # Pop 1
            [65, 1.3, 4e6]       # Pop 2
        ]
    elif Z==6:
        pars = [
            # a, gamma, R
            [100, 1.4, 1.2e5], # Pop 1
            [6, 1.3, 4e6],     # Pop 2
        ]
    elif Z==8:
        pars = [
            # a, gamma, R
            [130, 1.4, 1.2e5], # Pop 1
            [7, 1.3, 4e6],     # Pop 2
        ]
    elif Z==26:
        pars = [
            # a, gamma, R
            [60, 1.3, 1.2e5],       # Pop 1
            [2.3, 1.2, 4e6],        # Pop 2
            [0.025, 1.2, 1.3e9]     # Pop 3
        ]
    flux = np.zeros_like(E)
    for AGammaR_List in pars:
        flux += flux_single_population(AGammaR_List=AGammaR_List)
    return flux