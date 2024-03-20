import pandas as pd
import math
import glob
from utils import *
# import yaml

class CorsikaSamples:
    def __init__(self, path_list, josn_list, num_events_list, mc_power_index, energy_range, costh_range, sample_cylinder_surface, select_id, obs_level='det_level') -> None:
        self.path_list = path_list
        self.json_list = josn_list
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
        self.particles = self.load_pars_and_prim()
        self.SelectedParticles = self.pick_up_particles(select_id)
        self.primaries = self.get_prims()

    def load_pars_and_prim(self):
        particles = []
        shower_id_start = 0
        for ibatch in range(len(self.pars_paths)):
            cur_pars = pd.read_parquet(self.pars_paths[ibatch])
            prim = pd.read_json(self.json_list[ibatch])
            merge = cur_pars.merge(prim, on='shower', how='left')
            cur_pars['primary'] = merge['E']
            cur_pars.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            particles.append(cur_pars)
            # with open(self.prim_paths[ievent]) as f:
            #     cur_prim = yaml.safe_load(f)
            #     cur_prim = pd.DataFrame.from_dict(cur_prim)
        return pd.concat(particles)

    def get_prims(self):
        prims=[]
        shower_id_start = 0
        for ibatch in range(len(self.json_list)):
            cur_prims = pd.read_json(self.json_list[ibatch])
            cur_prims.shower *= shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            prims.append(cur_prims)
        return pd.concat(prims)

    def pick_up_particles(self, id):
         selected = self.particles[self.particles['pdg'].isin(id)]
         return selected

def weight(particles, weighter, emin, emax):
    particles['weight'] = particles['E'].apply(weighter,args=(emin,emax))
    return particles

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    plots['E'] = PlotContainer(xlabel='primary energy(GeV)', ylabel='Counts', logx=False, logy=False, figname=save_dir + 'primary_distribution_E.jpg', bins=np.linspace(100,100000,100))
    return plots

class GlobalSetting:
    path_list = glob.glob("/media/ineffablord/T7/siMu_atm/data/out_no_sea_obs_1e2to1e5GeV_E_-1/part?*/my_shower/")
    json_list = glob.glob("/media/ineffablord/T7/siMu_atm/data/out_no_sea_obs_1e2to1e5GeV_E_-1/part?*/Primaries.json")
    corsika_samples = CorsikaSamples(path_list=path_list, josn_list=json_list,
                                     num_events_list=[1000 for i in path_list], mc_power_index=-1,
                                     energy_range=[1e2, 1e5], costh_range=[0.9, 1], sample_cylinder_surface=math.pi*2000**2,
                                     select_id=[13,-13])
    save_dir = './save/'
    plots = setup_plots(save_dir)

if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.corsika_samples.particles)