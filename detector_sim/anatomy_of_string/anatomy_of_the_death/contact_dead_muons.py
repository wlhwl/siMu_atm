import pandas as pd
import numpy as np

if __name__ == '__main__':
    mu = pd.read_csv('../dead_event.csv', header=None)
    mu.columns = ['muonid']
    mu[['x','y','px','py','pz']] = np.nan
    mu.set_index('muonid', inplace=True)
    muonid = np.unique(mu.index)
    for i in muonid:
        muon_inject = pd.read_json('../../single_muon_jobs/job_' + str(i//10260) + '/mc_events.json')
        muon_inject = muon_inject[muon_inject['event_id'] == i]
        par_in_det = muon_inject['particles_at_detector'].to_dict()
        mu.loc[i, 'x'] = par_in_det[i%10260][0]['x']
        mu.loc[i, 'y'] = par_in_det[i%10260][0]['y']
        mu.loc[i, 'px'] = par_in_det[i%10260][0]['px']
        mu.loc[i, 'py'] = par_in_det[i%10260][0]['py']
        mu.loc[i, 'pz'] = par_in_det[i%10260][0]['pz']
        print(i)

    mu.to_csv('./whole_body.csv')