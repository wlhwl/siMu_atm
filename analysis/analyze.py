from config import *
from utils import *
import numpy as np

def Simu_spectrum_weight(e,emin,emax):
    nor = np.log(emax/emin)
    density_reciprocal = e
    return nor*density_reciprocal/(emax - emin)

def Weighter(e,emin,emax):
    #lne = np.log(e)
    spec_wei = Simu_spectrum_weight(e,emin,emax)
    return spec_wei

if __name__ == '__main__':
    myset = GlobalSetting()
    particles = myset.corsika_samples.primaries
    # weight
    particles = weight(particles,Weighter,myset.corsika_samples.energy_range[0],myset.corsika_samples.energy_range[1])
    print(particles)

    plots = myset.plots
    variables = list(plots.keys())
    for var in variables:
        pc = plots.get(var)
        to_draw = particles[var]
        bins = pc.bins if hasattr(pc, 'bins') else 30
        hist, bins, _ = pc.ax.hist(to_draw, bins, histtype="step", linewidth=2)
        pc.apply_settings()
        pc.savefig()
