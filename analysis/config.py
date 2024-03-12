import pandas as pd
import math
import glob
from utils import *
# import yaml

class CorsikaSamples:
    def __init__(self, path_list, num_events_list, mc_power_index, energy_range, costh_range, sample_cylinder_surface, obs_level='sea_level') -> None:
        self.path_list = path_list
        self.num_events_list = num_events_list
        self.mc_power_index = mc_power_index
        self.energy_range = energy_range
        self.costh_range = costh_range
        self.costh_range = costh_range
        self.sample_cylinder_surface = sample_cylinder_surface
        if (obs_level!='sea_level') and (obs_level!='det_level'):
            raise Exception('obs_level should be either "sea_level" or "det_level", not {}'.format(obs_level))
        self.obs_level = obs_level
        self.pars_paths = [p + 'particles_' + obs_level + '/particles.parquet' for p in self.path_list]
        self.prim_paths = [p + 'primary_det' for p in self.path_list]
        self.particles = self.load_pars_and_prim()

    def load_pars_and_prim(self):
        particles = []
        shower_id_start = 0
        for ibatch in range(len(self.pars_paths)):
            cur_pars = pd.read_parquet(self.pars_paths[ibatch])
            cur_pars.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            particles.append(cur_pars)
            # with open(self.prim_paths[ievent]) as f:
            #     cur_prim = yaml.safe_load(f)
            #     cur_prim = pd.DataFrame.from_dict(cur_prim)
        return pd.concat(particles)

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    plots['kinetic_energy'] = PlotContainer(xlabel=r'$E_{\mu}$ [GeV]', ylabel='Counts', logx=True, logy=False, figname=save_dir + 'muon_e_spectrum.pdf', bins=np.logspace(2, 4, 21))
    return plots


class GlobalSetting:
    path_list = glob.glob("/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/no_neu_bin/out_cut100/part?/my_shower/")
    corsika_samples = CorsikaSamples(path_list=path_list, num_events_list=[200 for i in path_list], mc_power_index=-2,\
             energy_range=[1e9, 1e10], costh_range=[0.9, 1], sample_cylinder_surface=math.pi*2000**2)
    
    save_dir = './save/'
    plots = setup_plots(save_dir)


if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.corsika_samples.particles)