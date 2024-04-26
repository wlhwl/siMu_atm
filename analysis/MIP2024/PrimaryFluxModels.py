import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from scipy.interpolate import interp1d
from utils import *

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
        

class PolyGonatoModel(PrimaryFluxModel):
    def __init__(self):
        super().__init__()
        self.yc, self.epc = -4.7, 1.87

        # Table of yz, phiz, Ez
        self.parameter_table = {
            # proton
            1: [-2.71, 8.73e-2, 4.5e6], 
            # He
            2: [-2.64, 5.71e-2, 9e6],
            # CNO, here take it as C
            6: [-2.68, 3.24e-2, 3.06e7],
            # Mg, here take it as O
            8: [-2.67, 3.16e-2, 6.48e7],
            # Fe
            26: [-2.58, 2.18e-2, 1.17e8]
        }
    
    def __call__(self, Z, E):
        # Energy unit: GeV
        yz, phiz, Ez = self.parameter_table[Z]
        yc, epc = self.yc, self.epc
        return phiz * (E/1e3)**yz * ( 1 + (E/Ez)**epc )**((yc-yz)/epc) / 1e3


if __name__=='__main__':
    plot_dir = './plots/'
    default_color_list = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    Z_list = [1, 2, 6, 8, 26]
    model_dict = {
        'GST3': GST3Model(),
        'GSF': GSFModel(),
        'poly-gonato': PolyGonatoModel(),
    }

    # Plot settings
    E_POWER_INDEX = 2.7
    x_values = np.logspace(3,8,51)
    model_total_flux = {}
    model_pc = {}
    model_color = {
        'GST3': default_color_list[0],
        'GSF': default_color_list[1],
        'poly-gonato': default_color_list[2],
    }
    for name, model in model_dict.items():
        PC = PlotContainer if name!= 'GSF' else RatioPlotContainer
        pc_cur = PC(logx=True, logy=True, xlabel=r'Primary Energy [GeV]', ylabel=rf'$E^{{{E_POWER_INDEX}}}dN/dE$ [$m^{{-2}}s^{{-1}}sr^{{-1}}GeV^{{{E_POWER_INDEX-1:.1f}}}$]', 
            figname=plot_dir+f"model_{name}.pdf", x_values = x_values, ylim=(1, 1e5)) 
            
        total_flux = np.zeros_like(x_values)
        for i, Z in enumerate(Z_list):
            flux = model(Z, x_values)
            # Draw flux for current Z
            pc_cur.ax.plot(x_values, flux * x_values**E_POWER_INDEX, label=name + f' Z={Z}', color=default_color_list[i+3], linestyle='--')
            total_flux += flux
        model_total_flux[name] = total_flux
        model_pc[name] = pc_cur

    for name_i, pc_cur in model_pc.items():
        for name_j, total_flux in model_total_flux.items():
            pc_cur.ax.plot(x_values, total_flux * x_values**E_POWER_INDEX, label=name_j+' Total', color=model_color[name_j], linewidth=2)
        
        model_pc['GSF'].insert_data(
            x_values=x_values, y_value=model_total_flux['GSF']*x_values**E_POWER_INDEX, ary_index=0, label=name_i
        )
        model_pc['GSF'].insert_data(
            x_values=x_values, y_value=model_total_flux[name_i]*x_values**E_POWER_INDEX, ary_index=1, label=name_i
        )
        if name_i!='GSF':
            # Save current model flux
            pc_cur.apply_settings(if_legend=True)
            pc_cur.savefig()
    
    pc = model_pc['GSF']
    pc.draw_ratio(f'Ratio to GSF', draw_error=False)
    pc.apply_settings(if_legend=True, ratio_ylim=(0.5, 1.5))
    pc.savefig()
