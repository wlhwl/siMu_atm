import numpy as np
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
            # nz_arr = df['nz'].to_numpy()
            # print(nz_arr)

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
    def __init__(self, path_list, josn_list, num_events_list, id, energy_range, costh_range, sample_cube_surface, select_id, obs_level='det_level') -> None:
        self.path_list = path_list
        self.json_list = josn_list
        self.num_events_list = num_events_list
        self.id = id
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
            cur_prims = pd.read_json(self.json_list[ibatch])
            cur_prims.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            prims.append(cur_prims)
        return pd.concat(prims)

    def pick_up_particles(self, id):
         selected = self.particles[self.particles['pdg'].isin(id)]
         return selected

def weight(particles, id, weighter, prim_num, energy_range, costh_range):
    emin = energy_range[0]
    emax = energy_range[1]
    cos_range = costh_range[1] - costh_range[0]
    particles['weight'] = particles['E'].apply(weighter, args=(id, emin,emax,prim_num)) * cos_range
    return particles

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    #plots['nz'] = PlotContainer(xlabel=r'$cos \theta$', ylabel=r'$dN/(dE\cdot dcos\theta)$', logx=False, logy=True, figname=save_dir + 'primaries_zenith_distribution.jpg', bins=np.linspace(0,1,50))
    #plots['E'] = PlotContainer(xlabel=r'$E_{primary}$', ylabel=r'$E^{2.6}dN/dE$', logx=True, logy=True, figname=save_dir + 'primary_all_distribution_2.6.jpg', bins=np.logspace(3,8,50))
    plots['kinetic_energy'] = PlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}]$', logx=True, logy=True, figname=save_dir + 'muon_spectrum.jpg', bins=np.logspace(2,8,50))
    return plots

class GlobalSetting:
    corsika_samples=[]
    prim_par = ['p','He','C','O','Fe']
    group_file = ['out_1-100T','out_100T-100P_lydmbadkyi','out_100T-100P_lydmsidklydmba','out_100T-100P_lydklydmsi']
    sim_group = [['out_1-100T',[1e3, 1e5],[0, 1],[[1000], [5000], [2000], [2000], [2000]]],
                 ['out_100T-100P_lydmbadkyi',[1e5, 1e8],[0.8, 1],[[2000],[2000],[2000],[2000],[2000]]],
                 ['out_100T-100P_lydmsidklydmba',[1e5,1e8],[0.4, 0.8],[[2000],[2000],[2000],[2000],[2000]]],
                 ['out_100T-100P_lydklydmsi',[1e5, 1e8],[0, 0.4],[[2000],[2000],[2000],[2000],[2000]]]]
   
    for i in [0,1,2,3,4]:
        for j in [0,1,2,3]:
            # full angle 1-100TeV
            cor_path_list = glob.glob('/media/ineffablord/T7/siMu_atm/data/' + prim_par[i] + '/' + sim_group[j][0] + '/part*/my_shower/')
            cor_json_list = glob.glob('/media/ineffablord/T7/siMu_atm/data/' + prim_par[i] + '/' + sim_group[j][0] + '/part*/Primaries.json')
            c_s0 = CorsikaSamples(path_list=cor_path_list,  josn_list=cor_json_list, 
                                num_events_list=np.multiply(np.ones_like(cor_path_list,np.double),sim_group[j][3][i]),
                                    id = i, energy_range=sim_group[j][1], costh_range=sim_group[j][2], sample_cube_surface=5000*16000,
                                        select_id=[13,-13])
            corsika_samples.append(c_s0)

    mp_path_list = glob.glob("/media/ineffablord/T7/siMu_atm/data/mupage/3.5-2.5/*.evt")
    mupage_samples = MupageSamples(path_list=mp_path_list, energy_range=(1, 1e5), costh_range=(0,math.cos(85)), extended_can_radius=300)

    save_dir = './picture/'
    plots = setup_plots(save_dir)

if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.corsika_samples.particles)
