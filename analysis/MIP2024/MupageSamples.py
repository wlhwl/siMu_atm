import numpy as np
import pandas as pd
import math

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