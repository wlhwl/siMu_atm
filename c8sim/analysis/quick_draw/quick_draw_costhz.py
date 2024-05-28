import sys
sys.path.append('./../')
from utils import *
import numpy as np
import pandas as pd


if __name__ == '__main__':
    save_dir = './../picture/'
    plots = RatioPlotContainer(xlabel=r'$\cos \theta$', ylabel=r'$\frac{dN}{d\Omega dS dt} [m^{-2}s^{-1}sr^{-1}]$', 
                          logx=False, logy=False, 
                          figname=save_dir + 'muon_zenith_dis.jpg', bins=np.linspace(0,1,50))
    
    cor = pd.read_parquet('/media/ineffablord/T7/siMu_atm/analysis/quick_draw/muon_cor.parquet')
    cor = cor.loc[cor.x**2 + cor.y**2 < 2e3**2]
    # cor = cor.loc[cor['z']==500]
    hist, bins = np.histogram(-cor['nz'], bins=plots.bins, weights=cor['weight']*1.02)#, density=True)
    hist = hist / np.diff(bins)

    mp = pd.read_parquet('/media/ineffablord/T7/siMu_atm/analysis/quick_draw/muon_mp.parquet')
    # mp = mp.loc[mp.x**2 + mp.y**2 < 2e3**2]
    mp = mp.loc[mp['z']==500]
    hist_mp, bins = np.histogram(-mp['nz'], bins=plots.bins, weights=mp['weight']/2/np.pi)#, density=True)
    hist_mp = hist_mp / np.diff(bins)
    
    plots.ax_main.plot((bins[:-1]+bins[1:])/2, hist, label='CORSIKA 8')
    plots.ax_main.plot((bins[:-1]+bins[1:])/2, hist_mp, label='MUPAGE')

    plots.insert_data((bins[:-1]+bins[1:])/2, hist, 1, 'ratio')
    plots.insert_data((bins[:-1]+bins[1:])/2, hist_mp, 0, 'ratio')
    plots.draw_ratio(draw_error=False)
    plots.ax_main.legend()
    plots.apply_settings(ratio_ylim=[0,5])
    plots.apply_settings
    plots.savefig()

    
    
    
