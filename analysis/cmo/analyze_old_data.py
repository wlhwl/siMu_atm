import numpy as np
import pandas as pd
import math
import glob
import re
import os
import json
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
    def __init__(self, path_list, num_events_list, id, energy_range, costh_range, primary_flux_model, primary_Z=1, 
             sample_radius=(285**2+2e3**2)**0.5, obs_level='det_level') -> None:
        self.primary_flux_model = primary_flux_model
        self.path_list = path_list
        self.num_events_list = num_events_list
        self.id = id
        self.energy_range = energy_range # GeV
        self.costh_range = costh_range
        self.primary_Z = primary_Z
        self.sample_radius = sample_radius # m
        self.obs_level = obs_level
        self.pars_paths = [p + '/particles/particles.parquet' for p in self.path_list]
        self.json_list = [p for p in self.path_list]
        self.particles, self.primaries = self.load_pars_and_prim()
        self.reweight_primary()
        self.muons = self.particles.loc[self.particles.pdg.abs()==13]

    def load_pars_and_prim(self):
        particles = []
        prims=[]
        shower_id_start = 0
        for ibatch in range(len(self.path_list)):
            try:
                with open(self.json_list[ibatch]) as f:
                    f = json.load(f)
                    cur_pars = pd.json_normalize(f, record_path='particles_at_detector', meta=['event_id'])
                    cur_pars['shower'] = cur_pars.event_id
                    cur_pars['pdg'] = cur_pars['pdgid']
                    cur_pars['energy'] = np.linalg.norm(cur_pars[['px','py','pz']], axis=1)
                    cur_prims = pd.json_normalize(f, record_path='particles_in', meta=['event_id'])
                    cur_prims['shower'] = cur_prims.event_id
                    cur_prims['E'] = np.linalg.norm(cur_prims[['px','py','pz']], axis=1)
                    cur_prims['pdg'] = cur_prims['pdgid']
            except:
                print(f'Error with file: {self.pars_paths[ibatch]}')
                continue
            merge = cur_pars.merge(cur_prims, on='shower', how='left')

            # Record primary information
            cur_pars['primary_energy'] = merge['E']
            cur_pars['primary_Z'] = self.primary_Z
            cur_prims.shower += shower_id_start
            cur_prims['Z'] = self.primary_Z
            cur_pars.shower += shower_id_start
            # shower_id_start += self.num_events_list[ibatch]
            shower_id_start += len(cur_prims) 
            particles.append(cur_pars)
            prims.append(cur_prims)
        return pd.concat(particles).set_index('shower'), pd.concat(prims).set_index('shower')
    
    def reweight_primary(self):
        # MC_spectrum_integral = np.log(self.energy_range[1]/self.energy_range[0])
        MC_spectrum_integral = (1 / self.energy_range[0] - 1 / self.energy_range[1])
        E_prim = self.primaries.E
        
        # MC primary energy spectrum PDF
        MC_energy_pdf = E_prim**-2 / MC_spectrum_integral 

        # MC primary PDF
        MC_pdf = MC_energy_pdf / (self.costh_range[1]-self.costh_range[0]) / (2*math.pi) / (math.pi * self.sample_radius**2)

        # reweight: target PDF divided by (MC PDF times num_events)
        # self.primaries['weight'] = 1 / MC_pdf / sum(self.num_events_list) * self.primary_flux_model(Z=self.primary_Z, E=E_prim)
        self.primaries['weight'] = 1 / MC_pdf / sum(self.num_events_list) * self.primary_flux_model(Z=self.primary_Z, E=E_prim)
        self.particles['weight'] = self.primaries["weight"].loc[self.particles.index]

class GlobalSetting:
    save_dir = './save-old-data/'
    os.makedirs(save_dir, exist_ok=True)

    # Load / Save Corsika8 muon
    muon_c8_path = save_dir + 'muon_C8.csv'
    if os.path.exists(muon_c8_path):
        muon_c8 = pd.read_csv(muon_c8_path).set_index('shower')
    else:
        # primary_flux_model = GST3Model()
        primary_flux_model = PolyGonatoModel()
        # primary_flux_model = GSFModel()
        prim_par = ['p','He','C','Mg','Fe']
        prim_Z = [1, 2, 6, 12, 26]
        group_file = ['out_1-100T','out_100T-100P_lydmbadkyi','out_100T-100P_lydmsidklydmba','out_100T-100P_lydklydmsi']
        sim_group = [['out_1-100T',[1e3, 1e5],[0, 1],[[1000], [5000], [2000], [2000], [2000]]],
                    ['out_100T-100P_lydmbadkyi',[1e5, 1e8],[0.8, 1],[[2000],[2000],[2000],[2000],[2000]]],
                    ['out_100T-100P_lydmsidklydmba',[1e5,1e8],[0.4, 0.8],[[2000],[2000],[2000],[2000],[2000]]],
                    ['out_100T-100P_lydklydmsi',[1e5, 1e8],[0, 0.4],[[2000],[2000],[2000],[2000],[2000]]]]
        
        sample_file_prefix = "/lustre/collider/mocen/project/hailing/data/atm_muon/dataStore/"
        sample_detail_list = [
            # file_name, Z, energy_range, num_events
            # protons
            ["sample_E1e3-1e4_02May2023/batch0/", 1, [1e3, 1e4], 2e8],
            ["sample_E1e4-1e5_01May2023/batch0/", 1, [1e4, 1e5], 2e6],
            ["sample_E1e5-1e6_12June2023/batch*/",1, [1e5, 1e6], 3e6],
            ["sample_E1e6-1e7_12June2023/batch*/",1, [1e6, 1e7], 3e5],
            # He
            ["sample_He_E1e3-1e4_19May2023/batch0/", 2, [1e3, 1e4], 1e8],
            ["sample_He_E1e4-1e5_19May2023/batch0/", 2, [1e4, 1e5], 5e6],
            ["sample_He_E1e5-1e6_19May2023/batch0/", 2, [1e5, 1e6], 2.5e6],
            ["sample_He_E1e6-1e8_19May2023/batch0/", 2, [1e6, 1e8], 4.98e5],
            # CNO
            ["sample_CNO_E1e4-5e4_2June2023/batch*/", 6, [1e4, 5e4], 8e7],
            ["sample_CNO_E5e4-1e5_2June2023/batch*/", 6, [5e4, 1e5], 5e6],
            ["sample_CNO_E1e5-1e6_2June2023/batch*/", 6, [1e5, 1e6], 1e6],
            ["sample_CNO_E1e6-1e7_20June2023/batch*/", 6, [1e6, 1e7], 1e5],
            ["sample_CNO_E1e7-1e8/", 6, [1e7, 1e8], 19900],
            # Mg
            ["sample_Mg_E3e4-1e5/batch*/", 12, [3e4, 1e5], 1e7],
            ["sample_Mg_E1e5-1e6/batch*/", 12, [1e5, 1e6], 1e6],
            ["sample_Mg_E1e6-1e7/batch*/", 12, [1e6, 1e7], 1e5],
            ["sample_Mg_E1e7-1e8/batch*/", 12, [1e7, 1e8], 11600],
            # Fe
            ["sample_Fe_E6e4-1e5/batch*/", 26, [6e4, 1e5], 2e6],
            ["sample_Fe_E1e5-1e6/batch*/", 26, [1e5, 1e6], 7.6e5],
            ["sample_Fe_E1e6-1e7/batch*/", 26, [1e6, 1e7], 9.2e4],
        ]
        corsika_samples=[]
        primary_c8 = []
        muon_c8 = []
        total_showers = 0
        for fname, Z, energy_range, num_events in sample_detail_list:
            print(fname)
            cor_path_list = glob.glob(sample_file_prefix + fname + "/mc_events*.json")
            cur_sample = CorsikaSamples(
                path_list=cor_path_list, primary_flux_model=primary_flux_model, num_events_list=[num_events], id=0, energy_range=energy_range, costh_range=[0,1], primary_Z=Z
            )
            corsika_samples.append(cur_sample)

            muons = cur_sample.muons
            muons.index += total_showers
            muon_c8.append(muons)
            primary = cur_sample.primaries
            primary.index += total_showers
            primary_c8.append(primary)
            total_showers += primary.index.max() + 1

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
        mp_path_list = glob.glob("/lustre/collider/mocen/project/hailing/MUPAGE/TRIDENT/evt/mupage-run_01.evt")
        mupage_samples = MupageSamples(path_list=mp_path_list, energy_range=(1, 1e5), costh_range=(0,math.cos(85/180*math.pi)), extended_can_radius=300)
        muon_mupage = mupage_samples.muons
        muon_mupage.to_csv(muon_mupage_path)

def draw_muon_spectrum(myset, muon_corsika, muon_mupage, costh_cut=0.95, bins=np.logspace(2,5,31)):
    muon_ls = [muon_corsika, muon_mupage]
    name_ls = ["CORSIKA", "MUPAGE"]
    pc = RatioPlotContainer(xlabel=r'$E_{\mu}$', ylabel=r'$EdN/dE [s^{-1}m^{-2}sr^{-1}]$', logx=True, logy=True, figname=myset.save_dir + 'separate_muon_spectrum.pdf')

    for i in range(2):
        muon, name = muon_ls[i], name_ls[i]
        color = default_color_list[i]

        # Select down-going muons
        costh_cut = 0
        muon = muon.loc[muon.nz<-1*costh_cut]
        ary, weights = muon.energy.to_numpy(), muon.weight.to_numpy()

        # Considering sr in weight
        # weights unit: GeV s-1 m-2 sr-1
        weights = weights / (2*math.pi*(1-costh_cut)) * ary
        bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

        # y_value: EdN/dE/dOmega/dS/dt
        y_value = bin_contents/np.diff(bins)
        hist, bins, _ = pc.ax.hist(bins[:-1], bins=bins, weights=y_value, label=name, histtype="step", linewidth=2, color=color)

        # Calculate errors
        y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
        y_err = y_err2**0.5 / np.diff(bins)
        pc.ax.errorbar((bins[1:]+bins[:-1])/2, y_value, yerr=y_err, fmt='none', color=color)
        pc.insert_data(
            x_values=(bins[:-1]+bins[1:])/2, 
            y_value=hist, ary_index=i, label='flux', y_err=y_err, color=color
        )

        # Draw primary composition for COSIKA sample
        if i==0:
            primary_list = muon.primary_Z.unique()
            primary_Z = muon.primary_Z.to_numpy()
            for j, Z in enumerate(primary_list):
                color = default_color_list[j+2]
                mask = (primary_Z==Z)
                tmp_ary, tmp_weights = ary[mask], weights[mask]

                # Draw hists
                bin_contents, bins = np.histogram(tmp_ary, bins=bins, weights=tmp_weights)

                # y_value: EdN/dE/dOmega/dS/dt
                x_values=(bins[:-1]+bins[1:])/2
                y_values = bin_contents/np.diff(bins)
                pc.ax.plot(x_values, y_values, drawstyle='steps-mid', color=color, label=name+f' Z={Z:.0f}', linestyle='--')
    # draw_depth_intensity_with_sealevel_intensity_and_energy_loss(pc.ax)
    # apply settings & draw ratio plot
    pc.draw_ratio(f'Mupage / CORSIKA', draw_error=True)
    pc.apply_settings(ratio_ylim=(0.1, 10), if_legend=False)
    pc.ax_ratio.set_yscale('log')
    pc.ax.legend(fontsize=7, loc=None)
    pc.savefig()


def restrict_muons(muon_list):
    for i in range(len(muon_list)):
        muon = muon_list[i]
        # Restrict muons to be in the upper surface
        muon = muon.loc[muon.z>499]
        # Restrict muon energy to be 100~1e5 GeV
        muon = muon.loc[(muon.energy>100) & (muon.energy<1e5)]
        muon_list[i] = muon
    return muon_list

if __name__ == '__main__':
    ##corsika
    myset = GlobalSetting()
    muon_corsika, muon_mupage = restrict_muons([myset.muon_c8, myset.muon_mupage])
    muon_corsika['nz'] = muon_corsika['pz'] / muon_corsika['energy']
   
    draw_muon_spectrum(myset=myset, muon_corsika=muon_corsika, muon_mupage=muon_mupage)