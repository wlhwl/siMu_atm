import math
from config import *
from utils import *
import numpy as np

if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    muon_corsika = myset.muon_c8
    muon_mupage = myset.muon_mupage
   

    # Draw vertical muon flux
    bins=np.logspace(2,5,30)
    pc = RatioPlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', logx=True, logy=True, figname=myset.save_dir + 'vertical_muon_spectrum.jpg', )

    # c8 downgoing muon
    costh_cut = 0.95
    muon = muon_corsika.loc[muon_corsika.nz<-1*costh_cut]
    ary, weights = muon.kinetic_energy.to_numpy(), muon.weight.to_numpy()
    # weights unit: GeV s-1 m-2 rad-1
    weights = weights / (2*math.pi*(1-costh_cut)) * ary
    bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)
    y_value = bin_contents/np.diff(bins)
    hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label="CORSICA", histtype="step", linewidth=2)
    y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
    y_err = y_err2**0.5 / np.diff(bins)
    pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none')
    pc.insert_data(
        x_values=(bins[:-1]+bins[1:])/2, 
        y_value=hist, ary_index=0, label='energy', y_err=y_err
    )

    # mupage downgoing muon
    muon = muon_mupage.loc[muon_mupage.nz<-1*costh_cut]
    ary, weights = muon.energy.to_numpy(), muon.weight.to_numpy()
    # weights unit: GeV s-1 m-2 rad-1
    weights = weights / (2*math.pi*(1-costh_cut)) * ary
    bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)
    y_value = bin_contents/np.diff(bins)
    hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label="Mupage", histtype="step", linewidth=2)
    y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
    y_err = y_err2**0.5 / np.diff(bins)
    pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none')
    pc.insert_data(
        x_values=(bins[:-1]+bins[1:])/2, 
        y_value=hist, ary_index=1, label='energy', y_err=y_err
    )

    pc.draw_ratio(f'Mupage / CORSIKA', draw_error=True)
    pc.apply_settings(ratio_ylim=(0., 5))
    legend_loc = pc.legend_loc if hasattr(pc, 'legend_loc') else None
    pc.ax.legend(fontsize=10, loc=legend_loc)
    pc.savefig()
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
