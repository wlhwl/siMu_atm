import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import uproot
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, default='0')
parser.add_argument('--output', type=str, default='0')
args = parser.parse_args()
filename = '../../ana_string/database/string_' + args.input + '.root'

hits = uproot.open(filename)["string_tree"]
hits = hits.arrays(library="pd")
hits['wavelength'] = 1240/hits['e0']
hits['EventId'] = hits['EventId'].apply(lambda x: x[0])
hits.set_index('EventId', inplace=True)

qe_data = pd.read_csv('../qe.csv')
qe = interp1d(qe_data['wavelength'], qe_data['qe'], kind='cubic', fill_value='extrapolate')

hits['qe'] = hits['wavelength'].apply(lambda x: np.random.random(len(x)) < qe(x)/100)

hits['t0'] = hits.apply(lambda row: [val for i, val in enumerate(row['t0']) if row['qe'][i]], axis=1)
hits['e0'] = hits.apply(lambda row: [val for i, val in enumerate(row['e0']) if row['qe'][i]], axis=1)
hits['wavelength'] = hits.apply(lambda row: [val for i, val in enumerate(row['wavelength']) if row['qe'][i]], axis=1)
hits['PmtId'] = hits.apply(lambda row: [val for i, val in enumerate(row['PmtId']) if row['qe'][i]], axis=1)
hits['DomId'] = hits.apply(lambda row: [val for i, val in enumerate(row['DomId']) if row['qe'][i]], axis=1)
hits = hits.drop('qe',axis=1)

hits = hits[hits['wavelength'].apply(len) > 0]

muon = pd.read_json('../muon.json')

hits.loc[:, 'weights'] = muon.loc[muon['event_id'].isin(hits.index), 'weights'].apply(lambda x: x.get('spectrum'))

# 创建一个新的列，表示weights值是否和上一行的值相同
hits.loc[:,'weights_shifted'] = hits['weights'].shift()
hits['weights_same_as_previous'] = hits['weights'] == hits['weights_shifted']

# 根据新的列对数据进行分组，并将每一组的数据合并到一起
hits['group'] = ((hits['weights_same_as_previous']==False)).cumsum()
merged_hits = hits
# # 删除不再需要的列
merged_hits = merged_hits.drop(['weights_shifted', 'weights_same_as_previous'], axis=1)

merged_hits.set_index('group', inplace=True)

time_window = 20
domhitname = 'domhits' + args.output + '.csv'
domhitfile = open(domhitname, 'a')
max_group = merged_hits.index.max()
for i in range(1, max_group+1):
    pmthits = np.zeros((21,31,20000))
    time_group = np.array(merged_hits.loc[i]['t0']).ravel().flatten()
    pmt_group = np.array(merged_hits.loc[i]['PmtId']).ravel().flatten()
    dom_group = np.array(merged_hits.loc[i]['DomId']).ravel().flatten()
    weights = np.array(merged_hits.loc[i]['weights']).ravel()[0]
    
  
    time_group = np.hstack(time_group)
    pmt_group = np.hstack(pmt_group)
    dom_group = np.hstack(dom_group)

    time_group = time_group - min(time_group)
    time_group = time_group//time_window

    dom_group = dom_group.astype(int)
    pmt_group = pmt_group.astype(int)
    time_group = time_group.astype(int)

    for j in range(len(time_group)):
        if ((time_group[j] < 20000) & (time_group[j] >= 0)):
            pmthits[dom_group[j]%21][pmt_group[j]][time_group[j]] = 1
            
    domhits = pmthits.sum(axis=1)
    for k in range(21):
        for l in range(1,21):
            count = np.count_nonzero(domhits[k] == l)
            domhitfile.write(str(count))
            domhitfile.write(',')
        domhitfile.write(str(weights))
        domhitfile.write('\n')
    
domhitfile.close()

    
