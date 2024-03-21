from config import *
from utils import *
import numpy as np

def Simu_spectrum_weight(e,emin,emax):
    nor = np.log(emax/emin)
    density_reciprocal = e
    return nor*density_reciprocal/(emax - emin)

def calculate_dnde(rc,a,gamma,z,e):
    dnde = a * np.power(e,-gamma-1) * np.exp(-e/(z*rc)) # dnde = phi / e
    return dnde

def Weighter(e,emin,emax):
    spec_wei = Simu_spectrum_weight(e,emin,emax)
    dnde = calculate_dnde(120*10**3,7000,1.661,1,e) +\
            calculate_dnde(4*10**6,150,1.4,1,e) +\
            calculate_dnde(1.3*10**9,1.4,1.4,1,e)
    return spec_wei*dnde#*np.power(e,2.6)

if __name__ == '__main__':
    myset = GlobalSetting()
    #get particles
    particles = myset.corsika_samples.SelectedParticles
    # weight
    particles = weight(particles,Weighter,myset.corsika_samples.energy_range[0],myset.corsika_samples.energy_range[1])
    print(particles)

    plots = myset.plots
    variables = list(plots.keys())
    for var in variables:
        pc = plots.get(var)
        to_draw = particles[var]
        bins = pc.bins if hasattr(pc, 'bins') else 30
        hist, bins, _ = pc.ax.hist(to_draw, bins, histtype="step", linewidth=2, weights=particles['weight'])
        pc.apply_settings()
        pc.savefig()
