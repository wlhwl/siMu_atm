import numpy as np
import pandas as pd
import math
import json
import os
from PrimaryFluxModels import *
from utils import *


class SeaLevelSample:
    # Constants
    RAW_POWER_INDEX = -1
    SAMPLE_AREA = math.pi * 500**2
    SAMPLE_SUFFIX = "mc_events.json"
    DETECT_RADIUS = 400 # m
    NOMINAL_MODEL_NAME="GSF"
    MODEL_DICT = {
        "PolyGonato": PolyGonatoModel(),
        "GST3": GST3Model(),
        "GSF": GSFModel(),
    }
    def __init__(self, primary_flux_model=PolyGonatoModel()) -> None:
        # Primary flux model
        self.primary_flux_model = primary_flux_model

        # sample paths
        # only protons are simulated
        path_prefix = "/lustre/collider/mocen/project/hailing/data/atm_muon/dataStore/sealevel/"
        self.sample_path_list = [
            # [energy range (GeV), num_events_per_batch, path_list]
            [[1e2, 1e3], 1e6, [path_prefix+f"E1e2-1e3_5June2023/separate/batch{i}/" for i in range(1)]], 
            [[1e3, 1e4], 1e6, [path_prefix+f"E1e3-1e4_5June2023/separate/batch{i}/" for i in range(1)]], 
            [[1e4, 1e5], 1e5, [path_prefix+f"E1e4-1e5_5June2023/separate/batch{i}/" for i in range(10)]], 
            [[1e5, 1e6], 1e4, [path_prefix+f"E1e5-1e6_5June2023/separate/batch{i}/" for i in range(100)]], 
            [[1e6, 1e7], 1e3, [path_prefix+f"E1e6-1e7_5June2023/batch{i}/" for i in range(50)]], 
            [[1e7, 1e8], 100, [path_prefix+f"E1e7-1e8_5June2023/batch{i}/" for i in range(20)]], 
        ]
        # load particles
        self.primaries, self.muons = self.load_samples()
    
    def load_samples(self):
        # Log for all samples
        all_prims, all_muons= [], []
        # global shower id of each batch
        shower_id_start = 0
        for energy_range, num_events_per_batch, sample_path in self.sample_path_list:
            print("Loading ", sample_path)
            # Log for current sample
            num_sample_events = 0
            sample_prims, sample_muons = [], []
            for batch in sample_path:
                # Load json file
                with open(batch+self.SAMPLE_SUFFIX) as f:
                    particles = json.load(f)
                primary = pd.json_normalize(particles, record_path='particles_in', meta=['event_id'])
                muons = pd.json_normalize(particles, record_path='particles_at_detector', meta=['event_id'])
                muons = muons.loc[muons['pdgid'].abs()==13]
                muons = muons.loc[muons.x**2+muons.y**2<self.DETECT_RADIUS**2]
                del particles
                # Get energy
                primary['energy'] = np.linalg.norm(primary[['px','py','pz']], axis=1)
                muons['energy'] = np.linalg.norm(muons[['px','py','pz']], axis=1)
                # Set index to be shower
                primary['shower'] = primary['event_id'] + shower_id_start
                muons['shower'] = muons['event_id'] + shower_id_start
                muons = muons.set_index('shower')
                primary = primary.set_index('shower').loc[muons.index.unique()] # only record primaries that contain informations
                # Increas shower_id_start and num_sample_events
                shower_id_start += num_events_per_batch
                num_sample_events += num_events_per_batch
                # Save current batch
                sample_prims.append(primary)
                sample_muons.append(muons)
                del primary, muons

            # Merge batch 
            sample_prims, sample_muons = pd.concat(sample_prims), pd.concat(sample_muons)
            # Reweight current sample
            sample_prims, sample_muons = self.reweight_particles(primaries=sample_prims, particles=sample_muons, energy_range=energy_range, num_sample_events=num_sample_events)
            sample_muons['weight'] = sample_prims["weight"].loc[sample_muons.index]
            # Save current sample
            all_prims.append(sample_prims)
            all_muons.append(sample_muons)
            del sample_prims, sample_muons
        all_prims, all_muons = pd.concat(all_prims), pd.concat(all_muons)
        return all_prims, all_muons

    def reweight_particles(self, primaries, particles, energy_range, num_sample_events):
        # primary energy
        E_prim = primaries['energy'].to_numpy()
        # Normalization factor for MC PDF
        MC_spectrum_integral = 0
        if self.RAW_POWER_INDEX==-1:
            MC_spectrum_integral = np.log(energy_range[1]/energy_range[0])
        else:
            MC_spectrum_integral = (self.RAW_POWER_INDEX+1) * (energy_range[1]**(self.RAW_POWER_INDEX+1) - energy_range[0]**(self.RAW_POWER_INDEX+1))
        
        # MC primary energy spectrum PDF
        MC_energy_pdf = E_prim**self.RAW_POWER_INDEX / MC_spectrum_integral 
        # MC primary PDF
        MC_pdf = MC_energy_pdf / (2*math.pi*1) / (self.SAMPLE_AREA)

        # truth_flux = np.zeros_like(E_prim)
        # weight unit: s-1 m-2
        for name, model in self.MODEL_DICT.items():
            # Get truth flux. Assume all particles are equal to proton
            truth_flux = np.zeros_like(E_prim)
            for Z in [1, 2, 6, 8, 26]:
                truth_flux += model(Z=Z, E=E_prim)
            primaries['weight_'+name] = 1 / MC_pdf / num_sample_events * truth_flux
            particles['weight_'+name] = primaries["weight_"+name].loc[particles.index]  / (math.pi * self.DETECT_RADIUS**2)
            primaries["weight_"+name] = primaries["weight_"+name] / (self.SAMPLE_AREA)
            if name==self.NOMINAL_MODEL_NAME:
                primaries['weight'] = primaries['weight_'+name]
                particles['weight'] = particles['weight_'+name]
        return primaries, particles
        
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
    def draw_vertical_muon_flux_sea_level(cls, muon, ax, color=None, label='CORSIKA8', E_power_index=2.7, costh_cut=0.88, bins=np.logspace(2,6,41), draw_aux_line=True):
        # Select muons
        muon = muon.loc[muon.pz.abs()/muon.energy>costh_cut]
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
            ax.plot(x_values, analytical_intensity*x_values**E_power_index, label="Analytical", linewidth=2)

            # Draw IceCube result
            ic_energy = np.logspace(3.8, 5.6, 19)
            ic_flux = cls.icecube_sea_level_vertical_muon_flux(ic_energy)
            ax.plot(ic_energy, ic_flux*ic_energy**E_power_index, label='IceCube Fit', linestyle='--', zorder=3, color='red', linewidth=2)
        return x_values, y_values, y_err


if __name__=='__main__':
    save_dir = './save/'
    plot_dir = './plots/'
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    # Load / Save muon and primary
    muon_path, primary_path = save_dir + 'sealevel_muon.csv', save_dir+'sealevel_primary.csv'
    if os.path.exists(muon_path) and os.path.exists(primary_path):
        muon = pd.read_csv(muon_path).set_index('shower')
        # primary = pd.read_csv(primary_path).set_index('shower')
    else:
        samples = SeaLevelSample()
        muon, primary = samples.muons, samples.primaries
        samples.muons.to_csv(muon_path)
        primary.to_csv(primary_path)

    # Draw muon flux

    # Compare different primary model
    # Plot settings
    E_power_index = 3.2
    pc = RatioPlotContainer( logx=True, logy=True, xlabel=r'$E_\mu$ [GeV]', ylabel=rf'$E_\mu^{{{E_power_index}}}dN/dE$ [$m^{{-2}}s^{{-1}}sr^{{-1}}GeV^{{{E_power_index-1:.1f}}}$]', 
            title='', figname=plot_dir+"sealevel_flux.pdf", )
    bins = np.logspace(2, 5.5, 31)
    x_values = (bins[1:]+bins[:1])/2

    # Relative value: analytical one
    analytical_intensity = SeaLevelSample.intensity_sea_level(x_values, zenith=0)
    relative_y_values = analytical_intensity*x_values**E_power_index
    relative_y_err = np.zeros_like(relative_y_values)

    for i, model_name in enumerate(list(SeaLevelSample.MODEL_DICT.keys())):
        color = default_color_list[i]

        # Insert relative value
        pc.insert_data(x_values=x_values, y_value=relative_y_values, ary_index=0, label=model_name, y_err=relative_y_err)

        # Get new weight
        muon['weight'] = muon["weight_"+model_name]
        # Draw flux
        x_values, y_values, y_err = SeaLevelSample.draw_vertical_muon_flux_sea_level(muon, ax=pc.ax_main, E_power_index=E_power_index, color=color, label=f'CORSIKA: {model_name}', bins=bins, draw_aux_line=(model_name==SeaLevelSample.NOMINAL_MODEL_NAME))
        pc.insert_data(
            x_values=x_values, y_value=y_values, ary_index=1, label=model_name, y_err=y_err, color=color
        )

    pc.draw_ratio(f'Ratio to Analytical', draw_error=True)
    pc.apply_settings(if_legend=True, ratio_ylim=(0.5, 2))
    pc.savefig()

            

