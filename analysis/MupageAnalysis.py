import pandas as pd
import math
import glob
from utils import *


class MupageSamples:
    def __init__(self, path_list, energy_range, costh_range, extended_can_radius) -> None:
        self.path_list = path_list
        self.energy_range = energy_range
        self.costh_range = costh_range
        self.extended_can_radius = extended_can_radius
        self.muons, self.livetime, self.livetime_error = self.load_muons()
        self.num_showers = len(self.muons.index.unique())
        self.define_plots_vars()
    
    def load_muons(self):
        muons = pd.DataFrame()
        livetime, livetime_error2 = 0., 0.
        total_shower = 0.
        for fpath in self.path_list:
            print('loading ', fpath)

            """
            get mc_events.json from an output of MUPAGE
            """
            with open(fpath, 'r') as file:
                first_row = file.readline()
                # get livetime & total n events
                second_row = file.readline().split(":")[1].split()
                livetime += float(second_row[0])
                livetime_error2 += float(second_row[1])**2
        
            # Read the file into a DataFrame
            df = pd.read_csv(fpath, 
                            delim_whitespace=True,  # Use whitespace as the delimiter
                            skiprows=2,  # Skip the first two lines
                            names=['showerId', 'multiplicity', 'muonId', 'x', 'y', 'z', 'nx', 'ny', 'nz', 'energy', 'relativeT', 'pdgId']) # t in s and relativeT in ns
            # ShowerId starts from 0
            df['showerId'] = df['showerId'] + total_shower - 1
            df['muonId'] -= 1
            # intrisic bugs in MUPAGE
            df.loc[df.multiplicity==1, 'nz'] *= -1
            total_shower += len(df.showerId.unique())

            # merge shower
            muons = pd.concat([muons, df])
        return muons.set_index('showerId'), livetime, livetime_error2**0.5

    def define_plots_vars(self):
        self.muons['R'] = np.linalg.norm(self.muons[['x', 'y']], axis=1)
        self.muons['costh'] = -1*self.muons['nz']
        
        self.shower_vars = pd.DataFrame(index=self.muons.index.unique())
        self.shower_vars['showerEnergy'] = self.muons.energy.groupby('showerId').sum()
        self.shower_vars['muonMultiplicity'] = self.muons.multiplicity.groupby('showerId').first()
        self.shower_vars['timeDuration'] = self.muons.relativeT.groupby('showerId').max() - self.muons.relativeT.groupby('showerId').min()
        


    

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    # Plots for muons
    energy_bin = np.logspace(0, 5, 51)
    plots['energy'] = PlotContainer(shower_vars=False, xlabel=r'$E_{\mu}$ [GeV]', ylabel=r'$EdN/dE\ [s^{-1}m^{-2}]$', logx=True, logy=True, figname=save_dir + 'muon_e_spectrum.pdf', bins=energy_bin, weights=1./(math.log(energy_bin[1]/energy_bin[0])))
    plots['R'] = PlotContainer(shower_vars=False, xlabel=r'$radius_{\mu}$ [m]', ylabel=r'$[s^{-1}m^{-2}]$', logx=False, logy=False, figname=save_dir + 'muon_radius.pdf', bins=np.linspace(0, 300, 31))
    plots['costh'] = PlotContainer(shower_vars=False, xlabel=r'$cos(\theta)$', ylabel=r'$[s^{-1}m^{-2}]$', logx=False, logy=False, figname=save_dir + 'muon_costh.pdf', bins=np.linspace(0, 1, 31))

    # Plots for showers
    energy_bin = np.logspace(0, 5, 51)
    plots['showerEnergy'] = PlotContainer(shower_vars=True, xlabel=r'$E_{shower}$ [GeV]', ylabel=r'$EdN/dE\ [s^{-1}m^{-2}]$', logx=True, logy=True, figname=save_dir + 'shower_e_spectrum.pdf', bins=energy_bin, weights=1./(math.log(energy_bin[1]/energy_bin[0])))
    plots['muonMultiplicity'] = PlotContainer(shower_vars=True, xlabel=r'Muon Multiplicity', ylabel=r'$[s^{-1}m^{-2}]$', logx=False, logy=True, figname=save_dir + 'shower_muon_multiplicity.pdf', bins=np.linspace(1, 20, 20))
    plots['timeDuration'] = PlotContainer(shower_vars=True, xlabel=r'Event Time Duration [ns]', ylabel=r'$[s^{-1}m^{-2}]$', logx=False, logy=True, figname=save_dir + 'shower_time_duration.pdf', bins=np.linspace(0, 1000, 31))
    return plots


class GlobalSetting:
    path_list = glob.glob("/lustre/collider/mocen/project/hailing/atmos_muon/generator_workspace/evt/*.evt")
    mupage_samples = MupageSamples(path_list=path_list, energy_range=(1, 1e5), costh_range=(math.cos(85), 1), extended_can_radius=300)
    
    save_dir = './save-mupage/'
    plots = setup_plots(save_dir)


if __name__ == '__main__':
    # basic settings
    myset = GlobalSetting()
    mupage_samples = myset.mupage_samples
    muons = mupage_samples.muons
    muons = muons.loc[muons.z == 400] # only consider upper surface muons
    # costh_cut = 0.95
    # muons = muons.loc[muons.costh >= costh_cut] # only consider down-going muons
    showers = mupage_samples.shower_vars

    # for plottings
    plots = myset.plots
    variables = list(plots.keys())

    # weights unit: [s-1m-2]
    print(f"Total livetime: {mupage_samples.livetime} s")
    weights = 1 / mupage_samples.livetime / (math.pi * mupage_samples.extended_can_radius**2) # / (1-costh_cut) / (2*math.pi)

    # loop over all variables
    for var in variables:
        pc = plots.get(var)
        to_draw = showers[var] if pc.shower_vars else muons[var]
        
        # weighted histogram
        bins = pc.bins if hasattr(pc, 'bins') else 30
        cur_weights = np.ones_like(to_draw) * (weights*pc.weights if hasattr(pc, 'weights') else 1)
        hist, _ = np.histogram(to_draw, bins=bins, weights=cur_weights)
        
        # statistical error
        counts, _ = np.histogram(to_draw, bins=bins,)
        yerr = np.zeros_like(counts)
        mask = counts>0
        yerr[mask] = 1 / np.sqrt(counts[mask]) * hist[mask]

        # draw
        x_values = (bins[1:] + bins[:-1]) / 2
        pc.ax.errorbar(x_values, hist, yerr=yerr, fmt='-o', linewidth=2)
    
        pc.apply_settings()
        pc.savefig()