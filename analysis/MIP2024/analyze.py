import math
from config import *
from utils import *
import numpy as np
import copy

def draw_muon_distribution(myset, muon):
    plots = {
        'energy': PlotContainer(logx=True, logy=True, xlabel=r'$E_\mu$ [GeV]', ylabel=rf'$E_\mu dN/dE$ [$m^{{-2}}s^{{-1}}sr^{{-1}}$]', 
            title='', figname=myset.plot_dir+"muon_spectrum.pdf", bins = np.logspace(2, 5, 41), var_name='energy', weight_scale_cz=1),
        'vertical energy spectrum': PlotContainer(logx=True, logy=True, xlabel=r'Vertical Muon $E_\mu$ [GeV]', ylabel=rf'$E_\mu dN/dE$ [$m^{{-2}}s^{{-1}}sr^{{-1}}$]', 
            title='', figname=myset.plot_dir+"muon_spectrum_vertical.pdf", bins = np.logspace(2, 5, 41), czcut=0.95, var_name='energy', weight_scale_cz=1/0.05),
        'nz': PlotContainer(logx=False, logy=True, xlabel=r'$cos(\theta)$', ylabel=rf'$dN/dcos\theta$ [$m^{{-2}}s^{{-1}}sr^{{-1}}$]', 
            title='', figname=myset.plot_dir+"muon_cosz.pdf", bins = np.linspace(-1, 0, 41), var_name='nz', weight_scale_cz=1) 
    }

    # Loop to draw
    for pc_name, pc in plots.items():
        bins = pc.bins
        if hasattr(pc, 'czcut'):
            cur_muon = muon.loc[muon['nz'].abs()>pc.czcut]
        else:
            cur_muon = muon

        primary_Z = cur_muon['primary_Z'].to_numpy()
        ary = cur_muon[pc.var_name].to_numpy()
        # weight unit: s-1 m-2 rad-1
        weight = (cur_muon['weight'].to_numpy()) * pc.weight_scale_cz / (2*math.pi) 
        # For *E*dN/dE
        if "energy" in pc_name:
            weight = weight * ary

        bin_contents, bins = np.histogram(ary, bins=bins, weights=weight)

        # Draw result
        color = default_color_list[0]
        y_values = bin_contents/np.diff(bins)
        x_values = (bins[1:]+bins[:-1])/2
        pc.ax.step(x_values, y_values, label="Total", where='mid', linewidth=2, color=color)
        # Calculate errors
        y_err2, bins = np.histogram(ary, bins=bins, weights=(weight)**2)
        y_err = y_err2**0.5 / np.diff(bins)
        pc.ax.errorbar(x_values, y_values, yerr=y_err, fmt='none', color=color)

        # Categorize into primary 
        for j, Z in enumerate(np.unique(primary_Z)):
            color = default_color_list[j+1]
            mask = (primary_Z==Z)
            tmp_ary, tmp_weight = ary[mask], weight[mask]
            bin_contents, bins = np.histogram(tmp_ary, bins=bins, weights=tmp_weight)
            x_values=(bins[:-1]+bins[1:])/2
            y_values = bin_contents/np.diff(bins)
            pc.ax.plot(x_values, y_values, drawstyle='steps-mid', color=color, label=f'Z={Z:.0f}', linestyle='--')
        pc.apply_settings(if_legend=True)
        pc.savefig()

def draw_bundle_var(myset, muon_corsika, muon_mupage):
    # Draw bundle variables
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]

    # Define vars to draw
    plots = {
        'multiplicity': RatioPlotContainer(xlabel='Multiplicity, m', ylabel=r'$dN/dm dS dt\ \ [s^{-1}m^{-2}]$', logx=False, logy=True, figname=myset.save_dir + 'bundle_multiplicity.pdf', bins=np.linspace(0.5, 40.5, 41), functor='count') ,
        'bund_energy': RatioPlotContainer(xlabel=r'Bundle Energy [GeV]', ylabel=r'$EdN/dE dS dt\ \ [s^{-1}m^{-2}]$', logx=True, logy=True, figname=myset.save_dir + 'bundle_energy.pdf', bins=np.logspace(2, 6, 41), functor='sum') ,
    }

    for var, pc in plots.items():
        bins = pc.bins
        functor = pc.functor
        for i, muon in enumerate(muon_ls):
            color = default_color_list[i]
            ary = muon.groupby('shower')['energy'].agg(functor).to_numpy()
            weights = muon.groupby('shower')['weight'].first().to_numpy()
            if 'energy' in var:
                weights = weights * ary

            # Insert array
            bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

            # y_value: dN/dx/dS/dt
            y_value = bin_contents/np.diff(bins)
            hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label=name_ls[i], histtype="step", linewidth=2, color=color)

            # Calculate errors
            y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
            y_err = y_err2**0.5 / np.diff(bins)
            pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none', color=color)
            pc.insert_data(
                x_values=(bins[:-1]+bins[1:])/2, 
                y_value=hist, ary_index=i, label='flux', y_err=y_err, color=color
            )

        # apply settings & draw ratio plot
        pc.draw_ratio(f'Mupage / CORSIKA', draw_error=True)
        pc.apply_settings(ratio_ylim=(0.1, 10), if_legend=False)
        pc.ax_ratio.set_yscale('log')
        legend_loc = None 
        pc.ax.legend(fontsize=7, loc=legend_loc)
        pc.savefig()

def restrict_muons(muon_list):
    for i in range(len(muon_list)):
        muon = muon_list[i]
        # Restrict muons to be in the upper surface
        muon = muon.loc[muon.z>499]
        # Restrict muon energy to be 100~1e5 GeV
        muon = muon.loc[(muon.energy>100) & (muon.energy<1e5)]
        muon_list[i] = muon
    return muon_list

if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    muon_corsika, muon_mupage = restrict_muons([myset.muon_c8, myset.muon_mupage])
   
    draw_muon_distribution(myset=myset, muon=muon_corsika)
    # draw_muon_spectrum(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # # Draw costh
    # draw_muon_costh(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # # Draw boundle variables
    # draw_bundle_var(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # Draw primary spectrum
    pc =  PlotContainer(xlabel=r'Primary Energy [GeV]', ylabel=r'$EdN/dE dS dt\ \ [s^{-1}m^{-2}]$', logx=True, logy=True, figname=myset.plot_dir + 'primary_spectrum.pdf') 
    bins=np.logspace(3, 7.5, 41)
    muon = muon_corsika
    ary = muon.groupby('shower')['primary_energy'].first().to_numpy()
    weights = muon.groupby('shower')['weight'].first().to_numpy() * ary 
    primary_Z = muon.groupby('shower')['primary_Z'].first().to_numpy()
    bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

    # y_value: dN/dx/dS/dt
    y_value = bin_contents/np.diff(bins)
    hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, histtype="step", linewidth=2, label='Total')

    # Calculate errors
    y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
    y_err = y_err2**0.5 / np.diff(bins)
    pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none')

    # Categorize into primary 
    for j, Z in enumerate(np.unique(primary_Z)):
        color = default_color_list[j+1]
        mask = (primary_Z==Z)
        tmp_ary, tmp_weights = ary[mask], weights[mask]
        bin_contents, bins = np.histogram(tmp_ary, bins=bins, weights=tmp_weights)
        x_values=(bins[:-1]+bins[1:])/2
        y_values = bin_contents/np.diff(bins)
        pc.ax.plot(x_values, y_values, drawstyle='steps-mid', color=color, label=f'Z={Z:.0f}', linestyle='--')
    pc.apply_settings(if_legend=True)
    pc.savefig()
