import pandas as pd
import numpy as np
from utils import *

if __name__ == '__main__':
    pc = PlotContainer(xlabel='hitted string num per muon', ylabel='counts', logx=False, logy=True, figname='hitted_string_at_least_3doms.jpg', bins=np.linspace(1,125,125))
    doms = pd.read_csv('../survive_event.csv', header=None)
    doms.set_index(0, inplace=True)
    index = np.unique(doms.index)
    counts=[]
    for i in index:
        # string_id = np.unique(doms.loc[i]//21) #no selection
        ## select string with at least n doms
        dom = np.unique(doms.loc[i])
        string_id = dom//21
        indices = np.where(np.bincount(string_id) > 2 )[0]
        string_id = np.unique(string_id[np.isin(string_id, indices)])
        counts.append(len(string_id))
    
    to_draw = np.histogram(counts, bins=np.linspace(1,125,125))
    pc.ax.hist(counts, bins=np.linspace(1,125,125), histtype='step', color='b', label='>= 3 doms per string')
    pc.ax.legend()
    pc.apply_settings()
    pc.savefig()

        
        
    

