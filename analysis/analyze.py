
from config import *
from utils import *

if __name__=='__main__':
    myset = GlobalSetting()
    plots = myset.plots
    variables = list(plots.keys())
    particles = myset.corsika_samples.particles

    for var in variables:
        pc = plots.get(var)
        to_draw = particles[var]
        bins = pc.bins if hasattr(pc, 'bins') else 30
        hist, bins, _ = pc.ax.hist(to_draw, bins, histtype="step", linewidth=2)
    
        pc.apply_settings()
        pc.savefig()

