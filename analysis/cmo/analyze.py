import math
from config import *
from utils import *
import numpy as np

def draw_muon_spectrum(myset, muon_corsika, muon_mupage, costh_cut=0.95, bins=np.logspace(2,5,31)):
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]
    pc = RatioPlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', logx=True, logy=True, figname=myset.save_dir + 'muon_spectrum.pdf')

    for i in range(2):
        muon, name = muon_ls[i], name_ls[i]
        # Restrict muons to be in the upper surface
        muon = muon.loc[muon.z>499]
        color = default_color_list[i]

        # Select down-going muons
        costh_cut = 0
        muon = muon.loc[muon.nz<-1*costh_cut]
        ary, weights = muon.energy.to_numpy(), muon.weight.to_numpy()

        # Considering sr in weight
        # weights unit: GeV s-1 m-2 rad-1
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
    legend_loc = pc.legend_loc if hasattr(pc, 'legend_loc') else None
    pc.ax.legend(fontsize=7, loc=legend_loc)
    pc.savefig()
    
if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    muon_corsika = myset.muon_c8
    muon_mupage = myset.muon_mupage
   
    draw_muon_spectrum(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)

    exit(0)
    plots = myset.plots
    variables = list(plots.keys())
    
    for var in variables:
        pc = plots.get(var)
        to_draw = muon_corsika[var]
        bins = pc.bins if hasattr(pc, 'bins') else 30
        
        hist_w, bin_w = np.histogram(to_draw, bins,weights=muon_corsika['weight']*muon_corsika['kinetic_energy']*2*np.pi)
        hist_w = hist_w / np.diff(bin_w)

        hist_raw, _ = np.histogram(to_draw, bins)
        err = hist_w / np.sqrt(hist_raw)

        hist_mp, bin_mp = np.histogram(mp_particles['energy'],bins,weights=mp_particles['energy']*mp_weight)
        hist_mp = hist_mp / np.diff(bin_mp)
        
        # pc.ax.plot((bin_w[1:]+bin_w[:-1])/2,hist_w,linestyle='--',label='sim')
        pc.ax.errorbar((bins[1:]+bins[:-1])/2, hist_w, yerr=err, label="corsika")
        pc.ax.plot((bins[1:]+bins[:-1])/2,hist_mp, label="mupage")
        pc.ax.legend()
        pc.apply_settings()
        pc.savefig()
