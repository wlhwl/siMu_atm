import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

if __name__ == '__main__':
    mu = pd.read_csv('./center_string.csv', header=None)
    mu = mu[mu.iloc[:,1] // 21==0]
    mu.set_index(0, inplace=True)
    mu['nz'] = np.nan
    mu['r'] = np.nan
    muonid = np.unique(mu.index)
    for i in muonid:
        muon_inject = pd.read_json('../../single_muon_jobs/job_' + str(i//10260) + '/mc_events.json')
        muon_inject = muon_inject[muon_inject['event_id'] == i]
        par_in_det = muon_inject['particles_at_detector'].to_dict()
        nz = par_in_det[i%10260][0]['pz'] / np.sqrt(par_in_det[i%10260][0]['px']**2 + par_in_det[i%10260][0]['py']**2 + par_in_det[i%10260][0]['pz']**2)
        r = np.sqrt(par_in_det[i%10260][0]['x']**2 + par_in_det[i%10260][0]['y']**2)
        mu.loc[i, 'nz'] = nz
        mu.loc[i, 'r'] = r

    mu.to_csv('./center_string1.csv')