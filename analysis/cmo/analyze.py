import math
from config import *
from utils import *
import numpy as np

def draw_muon_spectrum(myset, muon_corsika, muon_mupage, costh_cut=0.95, bins=np.logspace(2,5,31)):
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]
    pc = RatioPlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', logx=True, logy=True, figname=myset.save_dir + 'separate_muon_spectrum.pdf')

    for i in range(2):
        muon, name = muon_ls[i], name_ls[i]
        color = default_color_list[i]

        # Select down-going muons
        costh_cut = 0
        muon = muon.loc[muon.nz<-1*costh_cut]
        ary, weights = muon.energy.to_numpy(), muon.weight.to_numpy()

        # Considering sr in weight
        # weights unit: GeV s-1 m-2 sr-1
        weights = weights / (2*math.pi*(1-costh_cut)) * ary
        bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

        # y_value: EdN/dE/dOmega/dS/dt
        y_value = bin_contents/np.diff(bins)
        hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label=name, histtype="step", linewidth=2, color=color)

        # Calculate errors
        y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
        y_err = y_err2**0.5 / np.diff(bins)
        pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none', color=color)
        pc.insert_data(
            x_values=(bins[:-1]+bins[1:])/2, 
            y_value=hist, ary_index=i, label='flux', y_err=y_err, color=color
        )

        # Draw primary composition for COSIKA sample
        if i==0:
            primary_list = muon.primary_Z.unique()
            primary_Z = muon.primary_Z.to_numpy()
            for j, Z in enumerate(primary_list):
                color = default_color_list[j+2]
                mask = (primary_Z==Z)
                tmp_ary, tmp_weights = ary[mask], weights[mask]

                # Draw hists
                bin_contents, bins = np.histogram(tmp_ary, bins=bins, weights=tmp_weights)

                # y_value: EdN/dE/dOmega/dS/dt
                x_values=(bins[:-1]+bins[1:])/2
                y_values = bin_contents/np.diff(bins)
                pc.ax.plot(x_values, y_values, drawstyle='steps-mid', color=color, label=name+f' Z={Z:.0f}', linestyle='--')
    # apply settings & draw ratio plot
    pc.draw_ratio(f'Mupage / CORSIKA', draw_error=True)
    pc.apply_settings(ratio_ylim=(0.1, 10), if_legend=False)
    pc.ax_ratio.set_yscale('log')
    pc.ax.legend(fontsize=7, loc=None)
    pc.savefig()


def draw_muon_costh(myset, muon_corsika, muon_mupage, bins=np.linspace(0,1,41)):
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]

    pc = RatioPlotContainer(xlabel=r'cos$\theta$', ylabel=r'$dN/d\Omega dS dt\ \ [s^{-1}m^{-2}sr^{-1}]$', logx=False, logy=True, figname=myset.save_dir + 'separate_muon_costh.pdf')
    for i in range(2):
        muon, name = muon_ls[i], name_ls[i]
        color = default_color_list[i]
        ary, weights = -1*muon.nz.to_numpy(), muon.weight.to_numpy()

        # Considering sr in weight
        # weights unit: s-1 m-2 rad-1
        weights = weights / (2*math.pi)
        bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

        # y_value: dN/dcosth/dOmega/dS/dt
        y_value = bin_contents/np.diff(bins)
        hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label=name, histtype="step", linewidth=2, color=color)

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
    legend_loc = "upper left"
    pc.ax.legend(fontsize=7, loc=legend_loc)
    pc.savefig()
    

def draw_bundle_var(myset, muon_corsika, muon_mupage):
    # Draw bundle variables
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]

    # Define vars to draw
    plots = {
        'multiplicity': RatioPlotContainer(xlabel='Multiplicity, m', ylabel=r'$dN/dm dS dt\ \ [s^{-1}m^{-2}]$', logx=False, logy=True, figname=myset.save_dir + 'bundle_multiplicity.pdf', bins=np.linspace(1, 40, 40), functor='count') ,
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
   
    draw_muon_spectrum(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # Draw costh
    draw_muon_costh(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # Draw boundle variables
    draw_bundle_var(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    # Draw primary spectrum
    pc =  PlotContainer(xlabel=r'Primary Energy [GeV]', ylabel=r'$EdN/dE dS dt\ \ [s^{-1}m^{-2}]$', logx=True, logy=True, figname=myset.save_dir + 'C8_primary_spectrum.pdf') 
    bins=np.logspace(3, 7.5, 41)
    muon = muon_corsika
    ary = muon.groupby('shower')['primary_energy'].first().to_numpy()
    weights = muon.groupby('shower')['weight'].first().to_numpy() * ary 
    bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

    # y_value: dN/dx/dS/dt
    y_value = bin_contents/np.diff(bins)
    hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, histtype="step", linewidth=2)

    # Calculate errors
    y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
    y_err = y_err2**0.5 / np.diff(bins)
    pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none')
    pc.apply_settings()
    pc.savefig()
