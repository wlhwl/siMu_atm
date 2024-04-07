import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from scipy.interpolate import interp1d

class PrimaryFluxModel(ABC):
    def __init__(self) -> None:
        pass
    @abstractmethod
    def __call__(self, Z, E):
        ...

class GST3Model(PrimaryFluxModel):
    def __init__(self):
        super().__init__()
        # dN/dE = \Sum a*E^(-gamma-1) * exp(-E / Z / R)
        self.flux_single_population = lambda AGammaR_List, Z, E: AGammaR_List[0] * np.power(E,-AGammaR_List[1]-1) * np.exp(-E/(Z*AGammaR_List[2]))
        self.parameter_table = {
            # proton
            1: [
                # a, gamma, R
                [7e3, 1.661, 1.2e5], # Pop 1
                [150, 1.4, 4e6],     # Pop 2
                [1.4, 1.4, 1.3e9],   # Pop 3
            ],
            # He
            2: [
                # a, gamma, R
                [3200, 1.58, 1.2e5], # Pop 1
                [65, 1.3, 4e6]       # Pop 2
            ],
            # C
            6: [
                # a, gamma, R
                [100, 1.4, 1.2e5], # Pop 1
                [6, 1.3, 4e6],     # Pop 2
            ],
            # O
            8: [
                # a, gamma, R
                [130, 1.4, 1.2e5], # Pop 1
                [7, 1.3, 4e6],     # Pop 2
            ],
            # Fe
            26: [
                # a, gamma, R
                [60, 1.3, 1.2e5],       # Pop 1
                [2.3, 1.2, 4e6],        # Pop 2
                [0.025, 1.2, 1.3e9]     # Pop 3
            ]
        }
    
    def __call__(self, Z, E):
        parameter_list = self.parameter_table[Z]
        flux = np.zeros_like(E)
        for AGammaR_List in parameter_list:
            flux += self.flux_single_population(AGammaR_List=AGammaR_List, Z=Z, E=E)
        return flux


class GSFModel(PrimaryFluxModel):
    def __init__(self, file_path='./gsf_data_table.txt') -> None:
        super().__init__()
        # Load data from text file
        self.file_path = file_path
        data = pd.read_csv(file_path, delim_whitespace=True, skiprows=4, header=None)

        # Define column names based on your description
        column_names = ['x'] + [str(i) for i in range(1, 29)]
        data.columns = column_names

        # Get flux calculator with interpolation
        energy = data['x'].to_numpy()
        self.flux_calculator = {
            # Proton
            1: interp1d(energy, data['1'].to_numpy(), kind='linear', fill_value='extrapolate'),
            # He
            2: interp1d(energy, data['2'].to_numpy(), kind='linear', fill_value='extrapolate'),
            # C
            6: interp1d(energy, data[[str(i) for i in range(3, 7)]].sum(axis=1).to_numpy(), kind='linear', fill_value='extrapolate'),
            # O
            8: interp1d(energy, data[[str(i) for i in range(7, 11)]].sum(axis=1).to_numpy(), kind='linear', fill_value='extrapolate'),
            # Fe 
            26: interp1d(energy, data[[str(i) for i in range(11, 29)]].sum(axis=1).to_numpy(), kind='linear', fill_value='extrapolate'),
        }

    def __call__(self, Z, E):
        return self.flux_calculator[Z](E)
        


if __name__=='__main__':
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    energy = np.logspace(3, 8, 30)
    Z_list = [1, 2, 6, 8, 26]
    model_dict = {
        'GST3': GST3Model(),
        'GSF': GSFModel()
    }

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5), dpi=600, constrained_layout=True)
    for name, model in model_dict.items():
        linestyle = '--' if name=='GST3' else None
        total_flux = np.zeros_like(energy)
        for i, Z in enumerate(Z_list):
            # flux += model(Z, energy)
            flux = model(Z, energy)
            ax.plot(energy, flux * energy**2.6, label=name + f' Z={Z}', color=default_color_list[i], linestyle=linestyle)
            total_flux += flux
        ax.plot(energy, total_flux * energy**2.6, label=name+' Total', color=default_color_list[i+1], linewidth=2, linestyle=linestyle)
    ax.set_ylabel(r'Flux$\times\ E^{2.6}$')
    ax.set_xlabel(r'Energy [GeV]')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    fig.savefig('./save/model_comparison.pdf')
    
    
