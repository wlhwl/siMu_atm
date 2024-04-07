from utils import *
from PrimaryFluxModels import *
import math
import numpy as np
import pandas as pd


if __name__=='__main__':
    flux_energy_power_index = 2.6
    # flux_energy_power_index = 0

    # Load MC
    primary_MC_path = 'save/primary_C8.csv'
    primary_MC = pd.read_csv(primary_MC_path)
    primary_Z = primary_MC.Z.to_numpy()
    MC_energy, MC_weight = primary_MC.E.to_numpy(), primary_MC.weight.to_numpy()
    sample_radius2 = (285**2+2e3**2)
    MC_weight = MC_weight / (2*math.pi) / (math.pi * sample_radius2) * MC_energy**flux_energy_power_index


    # Setup flux model
    Z_list = [1, 2, 6, 8, 26]
    model = GSFModel()

    pc = RatioPlotContainer(xlabel=r'Energy [GeV]', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', logx=True, logy=True, figname='save/MC_Model_comparison.pdf')

    # x values
    bins = np.logspace(3, 8, 51)
    x_values = (bins[1:] + bins[:-1])/2

    total_y_model, total_y_MC = np.zeros_like(x_values), np.zeros_like(x_values)
    for i, Z in enumerate(Z_list):
        color = default_color_list[i]
        # Draw model
        flux = model(Z, x_values) * x_values**flux_energy_power_index
        pc.ax.plot(x_values, flux, drawstyle='steps-mid', color=color, label=f'Model Z={Z:.0f}', linestyle='--')
        total_y_model += flux

        # Draw MC
        mask = (primary_Z==Z)
        tmp_ary, tmp_weights = MC_energy[mask], MC_weight[mask]
        bin_contents, bins = np.histogram(tmp_ary, bins=bins, weights=tmp_weights)
        flux = bin_contents/np.diff(bins)
        pc.ax.plot(x_values, flux, drawstyle='steps-mid', color=color, label=f'MC Z={Z:.0f}', linestyle=None)
        total_y_MC += flux
    
    pc.ax.plot(x_values, total_y_model, drawstyle='steps-mid', color='black', linestyle='--', linewidth=1, label='Model total')
    pc.insert_data(
        x_values=x_values, y_value=total_y_model, ary_index=0, label='flux'
    )

    pc.ax.plot(x_values, total_y_MC, drawstyle='steps-mid', color='black', linestyle=None, linewidth=1, label='MC total')
    pc.insert_data(
        x_values=x_values, y_value=total_y_MC, ary_index=1, label='flux'
    )

    # apply settings & draw ratio plot
    pc.draw_ratio(f'MC/Model', draw_error=False)
    pc.apply_settings(ratio_ylim=(0.1, 10), if_legend=False)
    pc.ax_ratio.set_yscale('log')
    legend_loc = pc.legend_loc if hasattr(pc, 'legend_loc') else None
    pc.ax.legend(fontsize=7, loc=legend_loc)
    pc.savefig()
    
