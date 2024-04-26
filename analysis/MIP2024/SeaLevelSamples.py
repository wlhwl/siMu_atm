import numpy as np
import pandas as pd
import math
import glob
import re
import os
from PrimaryFluxModels import *

class CorsikaSealevelSamples:
    NOMINAL_MODEL_NAME="GSF"
    MODEL_DICT = {
        "GSF": GSFModel(),
        "GST3": GST3Model(),
        "PolyGonato": PolyGonatoModel()
    }
    def __init__(self, path_list, josn_list, num_events_list, id, energy_range, costh_range, primary_Z=1, 
             sample_radius=500, detect_radius=400, obs_level='det_level') -> None:
        self.path_list = path_list
        self.json_list = josn_list
        self.num_events_list = num_events_list
        self.id = id
        self.energy_range = energy_range # GeV
        self.costh_range = costh_range
        self.primary_Z = primary_Z
        self.sample_radius = sample_radius # m
        self.detect_radius = detect_radius
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
                cur_prims['E'] = cur_prims['E0']
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
        # weight unit: s-1 m-2
        for name, model in self.MODEL_DICT.items():
            self.primaries['weight_'+name] = 1 / MC_pdf / len(self.primaries) * model(Z=self.primary_Z, E=E_prim)
            # restrict particle radius
            self.particles = self.particles.loc[self.particles.x**2 + self.particles.y**2 < self.detect_radius**2]
            self.particles['weight_'+name] = self.primaries["weight_"+name].loc[self.particles.index] / (math.pi * self.detect_radius**2)
            self.primaries["weight_"+name] = self.primaries["weight_"+name] / (math.pi * self.sample_radius**2)
            if name==self.NOMINAL_MODEL_NAME:
                self.primaries['weight'] = self.primaries['weight_'+name]
                self.particles['weight'] = self.particles['weight_'+name]
        
    @staticmethod
    def intensity_sea_level(Emu, zenith: float):
        """
        Emu in GeV, zenith in rad
        return: dN/(dE dOmega dS dt) with unit 1/(m2 s sr GeV)
        """
        cs = math.cos(zenith)
        return 0.14e4 * Emu**-2.7 * ( 1/(1 + 1.1*Emu*cs/115) + 0.054/(1 + 1.1*Emu*cs/850))

    @staticmethod
    def icecube_sea_level_vertical_muon_flux(energy):
        # https://pos.sissa.it/301/317
        # fitting result from icecube with muon cosz>0/88 & muon log(E/GeV) in (3.8 5.6)
        # energy: GeV
        # return: s-1m-2sr-1GeV-1
        return 9e-13*(energy/50000)**(-3.74)

    @classmethod
    def draw_vertical_muon_flux_sea_level(cls, muon, ax, color=None, label='CORSIKA8', E_power_index=2.7, costh_cut=0.9, bins=np.logspace(2,6,41), draw_aux_line=True):
        # Select muons
        muon = muon.loc[muon['nz'].abs()>costh_cut]
        ary, weights = muon['energy'].to_numpy(), muon['weight'].to_numpy()

        # Considering sr in weight
        # weights unit: s-1 m-2 sr-1 GeV**E_power_index
        weights = weights / (2*math.pi*(1-costh_cut)) * ary**E_power_index
        bin_contents, bins = np.histogram(ary, bins=bins, weights=weights)

        # Draw result
        # y_value: EdN/dE/dOmega/dS/dt
        y_values = bin_contents/np.diff(bins)
        x_values = (bins[1:]+bins[:-1])/2
        ax.step(x_values, y_values, label=label, where='mid', linewidth=2, color=color)
        # Calculate errors
        y_err2, bins = np.histogram(ary, bins=bins, weights=(weights)**2)
        y_err = y_err2**0.5 / np.diff(bins)
        ax.errorbar(x_values, y_values, yerr=y_err, fmt='none', color=color)

        if draw_aux_line:
            # Draw Analytical function
            analytical_intensity = cls.intensity_sea_level(x_values, zenith=0)
            ax.plot(x_values, analytical_intensity*x_values**E_power_index, label="Analytical", linewidth=2, color=default_color_list[4])

            # Draw IceCube result
            ic_energy = np.logspace(3.8, 5.6, 19)
            ic_flux = cls.icecube_sea_level_vertical_muon_flux(ic_energy)
            ax.plot(ic_energy, ic_flux*ic_energy**E_power_index, label='IceCube Fit', linestyle='--', zorder=3, color='red', linewidth=2)
        return x_values, y_values, y_err

if __name__=='__main__':
    save_dir = './save/'
    plot_dir = './plots/'

    # Load / Save Corsika8 muon
    muon_c8_path = save_dir + 'C8_sealevel_muon.parquet'
    primary_path = save_dir + '/C8_sealevel_primary.parquet'

    if os.path.exists(muon_c8_path) and os.path.exists(primary_path):
        muons = pd.read_parquet(muon_c8_path).set_index('shower')
    else:
        prim_par = ['p','He','C','O','Fe']
        prim_Z = [1, 2, 6, 8, 26]
        sim_group = [
                    ['out_100G-100T_lydmjqdkyi',[1e2, 1e5],[0.9, 1],[[2000], [2000], [2000], [2000], [2000]]],
                    ['out_100T-100P_lydmjqdkyi',[1e5, 1e8],[0.9, 1],[[2000],[2000],[2000],[2000],[2000]]],
                ]
        corsika_samples=[]
        primary_c8 = []
        muons = []
        total_showers = 0
        for i, primary_name in enumerate(prim_par):
            for settings in sim_group:
                # full angle 1-100TeV
                cor_path_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/' + primary_name + '/' + settings[0] + '/part*/my_shower/')
                cor_json_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/bin/' + primary_name + '/' + settings[0] + '/part*/Primaries.json')
                c_s0 = CorsikaSealevelSamples(path_list=cor_path_list,  josn_list=cor_json_list,
                                    num_events_list=np.multiply(np.ones_like(cor_path_list,np.double),settings[3][i]),
                                        id = i, energy_range=settings[1], costh_range=settings[2], primary_Z=prim_Z[i])
                corsika_samples.append(c_s0)
                cur_muon = c_s0.muons
                cur_muon.index += total_showers
                muons.append(cur_muon)
                primary = c_s0.primaries
                primary.index += total_showers
                primary_c8.append(primary)
                total_showers += len(primary)
        muons = pd.concat(muons)
        primary_c8 = pd.concat(primary_c8)
        muons.to_parquet(muon_c8_path)
        # weight unit: [s-1 m-2 sr-1]
        primary_c8.to_parquet(primary_path)

    # Draw muon flux
    # Compare different primary model
    # Plot settings
    E_power_index = 3.2
    pc = RatioPlotContainer( logx=True, logy=True, xlabel=r'$E_\mu$ [GeV]', ylabel=rf'$E_\mu^{{{E_power_index}}}dN/dE$ [$m^{{-2}}s^{{-1}}sr^{{-1}}GeV^{{{E_power_index-1:.1f}}}$]', 
            title='', figname=plot_dir+"sealevel_flux.pdf", )
    bins = np.logspace(2, 5.5, 31)
    x_values = (bins[1:]+bins[:1])/2

    # Relative value: analytical one
    analytical_intensity = CorsikaSealevelSamples.intensity_sea_level(x_values, zenith=0)
    relative_y_values = analytical_intensity*x_values**E_power_index
    relative_y_err = np.zeros_like(relative_y_values)

    for i, model_name in enumerate(list(CorsikaSealevelSamples.MODEL_DICT.keys())):
        color = default_color_list[i]

        # Insert relative value
        pc.insert_data(x_values=x_values, y_value=relative_y_values, ary_index=0, label=model_name, y_err=relative_y_err)

        # Get new weight
        muons['weight'] = muons["weight_"+model_name] * 1.5
        # Draw flux
        x_values, y_values, y_err = CorsikaSealevelSamples.draw_vertical_muon_flux_sea_level(muons, ax=pc.ax_main, E_power_index=E_power_index, color=color, label=f'CORSIKA: {model_name}', bins=bins, draw_aux_line=(model_name=="PolyGonato"))
        pc.insert_data(
            x_values=x_values, y_value=y_values, ary_index=1, label=model_name, y_err=y_err, color=color
        )

    pc.draw_ratio(f'Ratio to Analytical', draw_error=True)
    pc.apply_settings(if_legend=True, ratio_ylim=(0., 2))
    pc.savefig()

            


