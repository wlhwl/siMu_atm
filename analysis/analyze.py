import math
from config import *
from utils import *
import numpy as np

para_table = [[[120*10**3,7000,1.661,1],[4*10**6,150,1.4,1],[1.3*10**9,1.4,1.4,1]],
                  [[120*10**3,3200,1.58,2],[4*10**6,65,1.3,2],[1.3*10**9,0,0,2]],
                  [[120*10**3,100,1.4,6],[4*10**6,6,1.3,6],[1.3*10**9,0,0,6]],
                  [[120*10**3,130,1.4,8],[4*10**6,7,1.3,8],[1.3*10**9,0,0,8]],
                  [[120*10**3,60,1.3,26],[4*10**6,2.3,1.2,26],[1.3*10**9,0.025,1.2,26]]]

def Simu_spectrum_weight(e,emin,emax,prim_num):
    nor = np.log(emax/emin)
    density_reciprocal = e
    return nor*density_reciprocal/prim_num

def calculate_dnde(e,rc,a,gamma,z):
    dnde = a * np.power(e,-gamma-1) * np.exp(-e/(z*rc)) # dnde = phi / e
    return dnde

def GST3_Weighter(e,id,emin,emax,prim_num):
    spec_wei = Simu_spectrum_weight(e,emin,emax,prim_num)
    dnde=0
    for i in range(0,3):
        dnde += calculate_dnde(e,para_table[id][i][0],para_table[id][i][1],para_table[id][i][2],para_table[id][i][3])
    return dnde*spec_wei   #*np.power(e,2.6)  *e*np.log(10)

if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    corsika = []
    for i in range(0,20):
        prim_num = np.sum(myset.corsika_samples[i].num_events_list)
        #get particles
        cor_particles = myset.corsika_samples[i].SelectedParticles # primaries 
        # cor_particles = cor_particles.loc[cor_particles.z == 500] 
        # cor_particles = cor_particles.loc[np.sqrt(np.square(cor_particles.x)+np.square(cor_particles.y)) < 2000]

        # weight
        cor_particles_w = weight(cor_particles, myset.corsika_samples[i].id, GST3_Weighter,
                                prim_num, myset.corsika_samples[i].energy_range,
                                myset.corsika_samples[i].costh_range)
        corsika.append(cor_particles_w)

    corsika = pd.concat(corsika)
    corsika.to_parquet('quick_draw/muon_cor.parquet')

    #mupage
    mp_particles = myset.mupage_samples.muons
    mp_weight = 1 / myset.mupage_samples.livetime / (math.pi * myset.mupage_samples.extended_can_radius**2)
    mp_particles['weight'] = mp_weight
    #print(mp_particles)
    mp_particles.to_parquet('quick_draw/muon_mp.parquet')
   
    # plots = myset.plots
    # variables = list(plots.keys())
    
    # for var in variables:
    #     pc = plots.get(var)
    #     to_draw = corsika[var]

    #     # energy = np.logspace(3, 8, 50)
    #     # dnde=0
    #     # for j in [0,1,2,3,4]:
    #     #     for i in range(0,3):
    #     #         dnde += calculate_dnde(energy,para_table[j][i][0],para_table[j][i][1],para_table[j][i][2],para_table[j][i][3])
    #     # pc.ax.plot(energy, dnde*np.power(energy,2.6), label="GST3 spectrum")

    #     bins = pc.bins if hasattr(pc, 'bins') else 30
        
    #     hist_w, bin_w = np.histogram(to_draw, bins,weights=corsika['weight']*corsika['kinetic_energy']*2*np.pi)
    #     hist_w = hist_w / np.diff(bin_w)

    #     hist_raw, _ = np.histogram(to_draw, bins)
    #     err = hist_w / np.sqrt(hist_raw)

    #     hist_mp, bin_mp = np.histogram(mp_particles['energy'],bins,weights=mp_particles['energy']*mp_weight)
    #     hist_mp = hist_mp / np.diff(bin_mp)
        
    #     # pc.ax.plot((bin_w[1:]+bin_w[:-1])/2,hist_w,linestyle='--',label='sim')
    #     pc.ax.errorbar((bins[1:]+bins[:-1])/2, hist_w, yerr=err, label="corsika")
    #     pc.ax.plot((bins[1:]+bins[:-1])/2,hist_mp, label="mupage")
    #     pc.ax.legend()
    #     pc.apply_settings()
    #     pc.savefig()
