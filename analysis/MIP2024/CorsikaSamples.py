import numpy as np
import pandas as pd
import math
import glob
import re
import os

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
