import math
from config import *
from utils import *
import numpy as np
from scipy import special

def Simu_spectrum_weight(e,emin,emax,prim_num):
    nor = np.log(emax/emin)
    density_reciprocal = e
    return nor*density_reciprocal/(emax-emin)/prim_num

def calculate_dnde(e,rc,a,gamma,z):
    dnde = a * np.power(e,-gamma-1) * np.exp(-e/(z*rc))# dnde = phi / e
    return dnde

def Weighter(e,emin,emax,prim_num):
    spec_wei = Simu_spectrum_weight(e,emin,emax,prim_num)
    para_table = [[120*10**3,7000,1.661,1],[4*10**6,150,1.4,1],[1.3*10**9,1.4,1.4,1]]
    dnde=0
    for i in range(0,3):
        dnde += calculate_dnde(e,para_table[i][0],para_table[i][1],para_table[i][2],para_table[i][3])/1.5295/np.power(10.,-5)
    return dnde*spec_wei*np.log(10)*e #/5000/5000 *np.power(e,2.6)  *e /area #*e/100000 #*np.power(e,2.6)

if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    prim_num = np.sum(myset.corsika_samples.num_events_list)
    #get particles
    cor_particles = myset.corsika_samples.primaries
    cor_particles['kinetic_energy'] = cor_particles.E
    # cor_particles = myset.corsika_samples.SelectedParticles
    # cor_particles = cor_particles.loc[cor_particles.z == 500]
    # weight
    cor_particles_w = weight(cor_particles,Weighter,prim_num,myset.corsika_samples.energy_range[0],myset.corsika_samples.energy_range[1])
    ##mupage
    mp_particles = myset.mupage_samples.muons
    mp_weight = 1 / myset.mupage_samples.livetime / (math.pi * myset.mupage_samples.extended_can_radius**2)#

    plots = myset.plots
    variables = list(plots.keys())
    for var in variables:
        pc = plots.get(var)
        to_draw = cor_particles_w[var]

        bins = pc.bins if hasattr(pc, 'bins') else 30

        hist_raw, bin_raw = np.histogram(to_draw, bins)

        hist_w, bin_w = np.histogram(to_draw, bins,weights=cor_particles_w['weight'])
        hist_w = hist_w / np.diff(bin_w)
        hist_mp, bin_mp = np.histogram(mp_particles['energy'],bins,weights=np.ones_like(mp_particles['energy'])*mp_weight)
        hist_mp = hist_mp / np.diff(bin_mp)
        # pc.ax.plot((bin_w[1:]+bin_w[:-1])/2,hist_w,linestyle='--',label='sim')
        pc.ax.errorbar((bins[1:]+bins[:-1])/2,hist_w,yerr=(hist_w/np.sqrt(hist_raw)), label="p 1-100TeV")
        # pc.ax.plot((bins[1:]+bins[:-1])/2,hist_mp)
        energy = np.logspace(3, 5, 20)
        para_table = [[120*10**3,7000,1.661,1],[4*10**6,150,1.4,1],[1.3*10**9,1.4,1.4,1]]
        dnde=0
        for i in range(0,3):
            dnde += calculate_dnde(energy,para_table[i][0],para_table[i][1],para_table[i][2],para_table[i][3])
        pc.ax.plot(energy, dnde * energy)
        pc.ax.legend()
        pc.apply_settings()
        pc.savefig()
