from utils import *
import numpy as np
import pandas as pd

if __name__ == '__main__':
    save_dir = './picture/'
    plots = RatioPlotContainer(xlabel=r'$E_{\mu} [GeV]$', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', 
                          logx=True, logy=True, 
                          figname=save_dir + 'muon_spectrum_com.jpg', bins=np.logspace(2,5,50))
    
    cor = pd.read_parquet('./muon_cor.parquet')
    mp = pd.read_parquet('./muon_mp.parquet')
    
    cor = cor.loc[cor['z']==500]
    cor = cor.loc[-cor['nz'] > 0.95]
    mp = mp.loc[mp['z']==400]
    mp = mp.loc[-mp['nz'] > 0.95]
    
    hist_cor, bin_cor = np.histogram(cor['kinetic_energy'],plots.bins,weights=cor['weight']*cor['kinetic_energy']/0.05)
    hist_cor = hist_cor / np.diff(bin_cor)
    hist_cor_raw, _ = np.histogram(cor['kinetic_energy'],plots.bins)
    err_cor = hist_cor / np.sqrt(hist_cor_raw)

    hist_mp, bin_mp = np.histogram(mp['energy'],plots.bins,weights=mp['weight']*mp['energy']/0.05/2/np.pi)
    hist_mp = hist_mp / np.diff(bin_mp)
    hist_mp_raw, _ = np.histogram(mp['energy'],plots.bins)
    err_mp = hist_mp / np.sqrt(hist_mp_raw)

    plots.ax_main.errorbar((bin_mp[:-1] + bin_mp[1:]) / 2,hist_mp,yerr=err_mp,label='MUPAGE')
    plots.ax_main.errorbar((bin_cor[:-1] + bin_cor[1:]) / 2,hist_cor,yerr=err_cor,label='CORSIKA')

    plots.insert_data((bin_cor[:-1] + bin_cor[1:]) / 2,hist_cor, 1,'ratio')
    plots.insert_data((bin_mp[:-1] + bin_mp[1:]) / 2,hist_mp, 0,'ratio')
    plots.draw_ratio(draw_error=False)
    
    plots.ax_main.legend()
    plots.apply_settings(ratio_ylim=[0,5])
    plots.savefig()

    
    
    