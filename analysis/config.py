import numpy as np
import pandas as pd
import math
import glob
from utils import *
# import yaml
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
            #df = df[df['energy'] > 200]
            #df = df[df['nz'] < -0.9]

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

class CorsikaSamples:
    def __init__(self, path_list, josn_list, num_events_list, energy_range, costh_range, sample_cube_surface, select_id, obs_level='det_level') -> None:
        self.path_list = path_list
        self.json_list = josn_list
        self.num_events_list = num_events_list
        self.energy_range = energy_range
        self.costh_range = costh_range
        self.sample_cube_surface = sample_cube_surface
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
            file_path = self.json_list[ibatch]
            with open(file_path, 'r') as file:
                file_content = file.read()
            file_content_fixed = file_content.replace('p+', '"p+"')
            prim =  pd.read_json(file_content_fixed)
            # prim = pd.read_json(self.json_list[ibatch])
            merge = cur_pars.merge(prim, on='shower', how='left')
            cur_pars['E'] = merge['E']
            cur_pars.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            particles.append(cur_pars)
        return pd.concat(particles)

    def get_prims(self):
        prims=[]
        shower_id_start = 0
        for ibatch in range(len(self.json_list)):
            file_path = self.json_list[ibatch]
            with open(file_path, 'r') as file:
                file_content = file.read()
            file_content_fixed = file_content.replace('p+', '"p+"')
            cur_prims =  pd.read_json(file_content_fixed)
            # cur_prims = pd.read_json(self.json_list[ibatch])
            cur_prims.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            prims.append(cur_prims)
        return pd.concat(prims)

    def pick_up_particles(self, id):
         selected = self.particles[self.particles['pdg'].isin(id)]
         return selected

def weight(particles, weighter, prim_num, emin, emax):
    particles['weight'] = particles['E'].apply(weighter, args=(emin,emax,prim_num))
    return particles

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    #plots['nz'] = PlotContainer(xlabel=r'$cos \theta$', ylabel=r'$dN/(dE\cdot dcos\theta)$', logx=False, logy=True, figname=save_dir + 'primaries_zenith_distribution.jpg', bins=np.linspace(0,1,50))
    plots['kinetic_energy'] = PlotContainer(xlabel=r'$E \mu$', ylabel=r'$EdN/dE$', logx=True, logy=True, figname=save_dir + 'muon_e_distribution.jpg', bins=np.logspace(2,5,30))
    return plots

class GlobalSetting:
    cor_path_list = glob.glob("/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/p/out_1-100T/part?*/my_shower/")
    cor_json_list = glob.glob("/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/p/out_1-100T/part?*/Primaries.json")
    corsika_samples = CorsikaSamples(path_list=cor_path_list,  josn_list=cor_json_list,
                                     num_events_list=[5000 for i in cor_path_list],
                                     energy_range=[1e3, 1e5], costh_range=[0, 1], sample_cube_surface=5000*16000,
                                     select_id=[13,-13])
    mp_path_list = glob.glob("/lustre/collider/mocen/project/hailing/atmos_muon/generator_workspace/evt/*.evt")
    mupage_samples = MupageSamples(path_list=mp_path_list, energy_range=(1, 1e5), costh_range=(0,math.cos(85)), extended_can_radius=300)

    save_dir = './test_weight/'
    plots = setup_plots(save_dir)

if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.corsika_samples.particles)
