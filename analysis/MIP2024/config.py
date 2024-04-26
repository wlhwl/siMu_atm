import numpy as np
import pandas as pd
import math
import glob
import re
import os
from utils import *
from PrimaryFluxModels import *
from MupageSamples import MupageSamples
from CorsikaSamples import CorsikaSamples


class GlobalSetting:
    save_dir = './save/'
    plot_dir = './plots/'
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    # Load / Save Corsika8 muon
    muon_c8_path = save_dir + 'C8_muon.csv'
    if os.path.exists(muon_c8_path):
        muon_c8 = pd.read_csv(muon_c8_path).set_index('shower')
    else:
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
                cor_path_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/underseabin/' + primary_name + '/' + settings[0] + '/part*/my_shower/')
                cor_json_list = glob.glob('/lustre/neutrino/huangweilun/atmos_muon/COR_atm_muon/test_proton/underseabin/' + primary_name + '/' + settings[0] + '/part*/Primaries.json')
                c_s0 = CorsikaSamples(path_list=cor_path_list,  josn_list=cor_json_list,
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
        muon_c8.to_csv(muon_c8_path)
        # weight unit: [s-1 m-2 sr-1]
        primary_c8.to_csv(save_dir + '/C8_primary.csv')


    # Load / Save mupage muon
    muon_mupage_path = save_dir + 'mupage_muon.csv'
    if os.path.exists(muon_mupage_path):
        muon_mupage = pd.read_csv(muon_mupage_path).set_index('shower')
    else:
        mp_path_list = glob.glob("/lustre/collider/mocen/project/hailing/atmos_muon/generator_workspace/evt/*.evt")
        mupage_samples = MupageSamples(path_list=mp_path_list, energy_range=(1, 1e5), costh_range=(0,math.cos(85/180*math.pi)), extended_can_radius=300)
        muon_mupage = mupage_samples.muons
        muon_mupage.to_csv(muon_mupage_path)

if __name__=='__main__':
    myset = GlobalSetting()
    print(myset.muon_c8)
    print(myset.muon_mupage)
