import numpy as np
import pandas as pd
import math
import glob
import re
import os
from utils import *
from PrimaryFluxModels import *

class MupageSamples:
    def __init__(self, path_list, energy_range, costh_range, extended_can_radius) -> None:
        self.path_list = path_list
        self.energy_range = energy_range
        self.costh_range = costh_range
        self.extended_can_radius = extended_can_radius
        self.muons, self.livetime, self.livetime_error = self.load_muons()
        self.reweight()
        self.num_showers = len(self.muons.index.unique())
        self.define_plots_vars()

    def reweight(self):
        # weight unit: [s-1 m-2]
        weight = 1 / self.livetime / (math.pi * self.extended_can_radius**2)
        # restrict muons within certain region
        self.muons = self.muons.loc[(self.muons.x**2+self.muons.y**2<(self.extended_can_radius**2))]
        self.muons['weight'] = weight

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
                             names=['shower', 'multiplicity', 'muonId', 'x', 'y', 'z', 'nx', 'ny', 'nz', 'energy', 'relativeT', 'pdgId']) # t in s and relativeT in ns
            # remove rows with NaN values
            df.dropna(inplace=True)
            df['nz'] = df['nz'].astype(float)
            # ShowerId starts from 0
            df['shower'] = df['shower'] + total_shower - 1
            df['muonId'] -= 1
            # intrisic bugs in MUPAGE
            df.loc[df.multiplicity==1, 'nz'] *= -1

            total_shower += len(df.shower.unique())
            # merge shower
            muons = pd.concat([muons, df])
        return muons.set_index('shower'), livetime, livetime_error2**0.5

    def define_plots_vars(self):
        self.muons['R'] = np.linalg.norm(self.muons[['x', 'y']], axis=1)
        self.muons['costh'] = -1*self.muons['nz']

        self.shower_vars = pd.DataFrame(index=self.muons.index.unique())
        self.shower_vars['showerEnergy'] = self.muons.energy.groupby('shower').sum()
        self.shower_vars['muonMultiplicity'] = self.muons.multiplicity.groupby('shower').first()
        self.shower_vars['timeDuration'] = self.muons.relativeT.groupby('shower').max() - self.muons.relativeT.groupby('shower').min()

class CorsikaSamples:
    def __init__(self, path_list, josn_list, num_events_list, id, energy_range, costh_range, primary_flux_model, primary_Z=1, 
             sample_radius=(285**2+2e3**2)**0.5, obs_level='det_level') -> None:
        self.primary_flux_model = primary_flux_model
        self.path_list = path_list
        self.json_list = josn_list
        self.num_events_list = num_events_list
        self.id = id
        self.energy_range = energy_range # GeV
        self.costh_range = costh_range
        self.primary_Z = primary_Z
        self.sample_radius = sample_radius # m
        self.obs_level = obs_level
        self.pars_paths = [p + 'particles_' + obs_level + '/particles.parquet' for p in self.path_list]
        self.particles, self.primaries = self.load_pars_and_prim()
        self.reweight_primary()
        self.muons = self.particles.loc[self.particles.pdg.abs()==13]

    def load_pars_and_prim(self):
        particles = []
        prims=[]
        shower_id_start = 0
        for ibatch in range(len(self.pars_paths)):
            try:
                cur_pars = pd.read_parquet(self.pars_paths[ibatch])
                file_path = self.json_list[ibatch]
                with open(file_path, 'r') as file:
                    file_content = file.read()
                file_content_fixed = re.sub(r'"pdg":.*\n', '', file_content)
                cur_prims =  pd.read_json(file_content_fixed)
            except:
                print(f'Error with file: {self.pars_paths[ibatch]}')
                continue
            merge = cur_pars.merge(cur_prims, on='shower', how='left')

            cur_pars['energy'] = cur_pars.kinetic_energy + 0.1
            # Record primary information
            cur_pars['primary_energy'] = merge['E']
            cur_pars['primary_Z'] = self.primary_Z
            cur_prims.shower += shower_id_start
            cur_prims['Z'] = self.primary_Z
            cur_pars.shower += shower_id_start
            shower_id_start += self.num_events_list[ibatch]
            particles.append(cur_pars)
            prims.append(cur_prims)
        return pd.concat(particles).set_index('shower'), pd.concat(prims).set_index('shower')
    
    def reweight_primary(self):
        MC_spectrum_integral = np.log(self.energy_range[1]/self.energy_range[0])
        E_prim = self.primaries.E
        
        # MC primary energy spectrum PDF
        MC_energy_pdf = E_prim**-1 / MC_spectrum_integral 

        # MC primary PDF
        MC_pdf = MC_energy_pdf / (self.costh_range[1]-self.costh_range[0]) / (2*math.pi) / (math.pi * self.sample_radius**2)

        # reweight: target PDF divided by (MC PDF times num_events)
        # self.primaries['weight'] = 1 / MC_pdf / sum(self.num_events_list) * self.primary_flux_model(Z=self.primary_Z, E=E_prim)
        self.primaries['weight'] = 1 / MC_pdf / len(self.primaries) * self.primary_flux_model(Z=self.primary_Z, E=E_prim)
        self.particles['weight'] = self.primaries["weight"].loc[self.particles.index]

def setup_plots(save_dir:str, plots:dict=None):
    if plots==None:
        plots = {}
    #plots['nz'] = PlotContainer(xlabel=r'$cos \theta$', ylabel=r'$dN/(dE\cdot dcos\theta)$', logx=False, logy=True, figname=save_dir + 'primaries_zenith_distribution.jpg', bins=np.linspace(0,1,50))
    #plots['E'] = PlotContainer(xlabel=r'$E_{primary}$', ylabel=r'$E^{2.6}dN/dE$', logx=True, logy=True, figname=save_dir + 'primary_all_distribution_2.6.jpg', bins=np.logspace(3,8,50))
    plots['kinetic_energy'] = PlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}]$', logx=True, logy=True, figname=save_dir + 'muon_spectrum.jpg', bins=np.logspace(2,8,50))
    return plots

class GlobalSetting:
    save_dir = './save/'
    os.makedirs(save_dir, exist_ok=True)

    # Load / Save Corsika8 muon
    muon_c8_path = save_dir + 'muon_C8.csv'
    if os.path.exists(muon_c8_path):
        muon_c8 = pd.read_csv(muon_c8_path).set_index('shower')
    else:
        # primary_flux_model = GST3Model()
        # primary_flux_model = PolyGonatoModel()
        primary_flux_model = GSFModel()
        prim_par = ['p','He','C','O','Fe']
        prim_Z = [1, 2, 6, 8, 26]
        group_file = ['out_1-100T','out_100T-100P_lydmbadkyi','out_100T-100P_lydmsidklydmba','out_100T-100P_lydklydmsi']
        sim_group = [['out_1-100T',[1e3, 1e5],[0, 1],[[1000], [5000], [2000], [2000], [2000]]],
                    ['out_100T-100P_lydmbadkyi',[1e5, 1e8],[0.8, 1],[[2000],[2000],[2000],[2000],[2000]]],
                    ['out_100T-100P_lydmsidklydmba',[1e5,1e8],[0.4, 0.8],[[2000],[2000],[2000],[2000],[2000]]],
                    ['out_100T-100P_lydklydmsi',[1e5, 1e8],[0, 0.4],[[2000],[2000],[2000],[2000],[2000]]]]
        corsika_samples=[]
        primary_c8 = []
        muon_c8 = []
        total_showers = 0
        for i, primary_name in enumerate(prim_par):
            for settings in sim_group:
                # full angle 1-100TeV
                cor_path_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/' + primary_name + '/' + settings[0] + '/part*/my_shower/')
                cor_json_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/' + primary_name + '/' + settings[0] + '/part*/Primaries.json')
                c_s0 = CorsikaSamples(path_list=cor_path_list,  josn_list=cor_json_list, primary_flux_model=primary_flux_model,
                                    num_events_list=np.multiply(np.ones_like(cor_path_list,np.double),settings[3][i]),
                                        id = i, energy_range=settings[1], costh_range=settings[2], primary_Z=prim_Z[i])
                corsika_samples.append(c_s0)
                muons = c_s0.muons
                muons.index += total_showers
                muon_c8.append(muons)
                primary = c_s0.primaries
                primary.index += total_showers
                primary_c8.append(primary)
                total_showers += len(primary)
        muon_c8 = pd.concat(muon_c8)
        primary_c8 = pd.concat(primary_c8)
        # restrict muons within certain region
        muon_c8 = muon_c8.loc[(muon_c8.x**2+muon_c8.y**2<(2e3**2))]
        # weight unit: [s-1 m-2]
        muon_c8['weight'] = muon_c8.weight / (math.pi * ( 2e3**2) )

        muon_c8.to_csv(muon_c8_path)
        # weight unit: [s-1 m-2 sr-1]
        primary_c8.to_csv(save_dir + '/primary_C8.csv')


    # Load / Save mupage muon
    muon_mupage_path = save_dir + 'muon_mupage.csv'
    if os.path.exists(muon_mupage_path):
        muon_mupage = pd.read_csv(muon_mupage_path).set_index('shower')
    else:
        mp_path_list = glob.glob("/lustre/collider/mocen/project/hailing/atmos_muon/generator_workspace/evt/*.evt")
        mupage_samples = MupageSamples(path_list=mp_path_list, energy_range=(1, 1e5), costh_range=(0,math.cos(85/180*math.pi)), extended_can_radius=300)
        muon_mupage = mupage_samples.muons
        muon_mupage.to_csv(muon_mupage_path)
    plots = setup_plots(save_dir)

if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.muon_c8)
    print(myset.muon_mupage)
