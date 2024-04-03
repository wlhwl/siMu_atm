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
        dnde += calculate_dnde(e,para_table[i][0],para_table[i][1],para_table[i][2],para_table[i][3])/1.808/np.power(10.,-5)
    return dnde*spec_wei*np.log(10)#*e*np.power(e,2.6)  *e /area #*e/100000 #*np.power(e,2.6)

def e_1(e):
    value = np.power(e,-1)/np.log(10**5/10**2)*100000
    return value

if __name__ == '__main__':
    myset = GlobalSetting()
    prim_num = np.sum(myset.corsika_samples.num_events_list)
    prim_num_add = np.sum(myset.corsika_samples_add.num_events_list)
    #get particles
    particles = myset.corsika_samples.SelectedParticles
    particles_add = myset.corsika_samples_add.SelectedParticles
    # weight
    particles_w = weight(particles,Weighter,prim_num,myset.corsika_samples.energy_range[0],myset.corsika_samples.energy_range[1])
    particles_w_add = weight(particles_add,Weighter,prim_num_add,myset.corsika_samples_add.energy_range[0],myset.corsika_samples_add.energy_range[1])

    plots = myset.plots
    variables = list(plots.keys())
    for var in variables:
        pc = plots.get(var)
        to_draw = particles_w[var]
        to_draw_add = particles_w_add[var]

        bins = pc.bins if hasattr(pc, 'bins') else 30

        # rawhist, bin_r = np.histogram(to_draw, bins)
        # rawhist_add, bin_s = np.histogram(to_draw_add, bins)
        # rawhist = rawhist*0.1 + rawhist_add*0.9
        # rawhist = rawhist / np.diff(bin_r)
        # errhist, bin_s = np.histogram(to_draw, bins, weights=particles_w['weight'])
        # errhist_add, bin_r = np.histogram(to_draw_add,bins,weights=particles_w_add['weight'])
        # errhist = errhist + errhist_add
        # errhist = errhist / np.diff(bin_s)
        hist_w, bins = np.histogram(to_draw, bins, weights=particles_w['weight'])
        hist_w = hist_w / np.diff(bins) * 0.1
        hist_w_add, bin_add = np.histogram(-to_draw_add,bins, weights=particles_w_add['weight'])
        hist_w_add = hist_w_add / np.diff(bin_add) * 0.9
        hist = hist_w_add + hist_w

        pc.ax.plot((bins[1:]+bins[:-1])/2,hist,linestyle='--',label='sim')
        #pc.ax.errorbar((bins[1:]+bins[:-1])/2,hist,yerr=(errhist/np.sqrt(rawhist)))
        #pc.ax.plot(bins,e_1(bins) ,linewidth=2,label='math')
        pc.ax.legend()
        pc.apply_settings()
        pc.savefig()
