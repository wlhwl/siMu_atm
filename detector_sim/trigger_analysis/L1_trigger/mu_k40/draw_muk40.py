from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from time_estimate import *
from utils import *
import numpy as np
import pandas as pd

def k40(pc_cl):
    num_per_run = 20000000
    k40coin_num = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    save_dir = './'
    world_radius = 30
    csv = pd.read_csv(save_dir + '30m.csv')
    event_num = num_per_run * len(csv)
    time = time_estimate(event_num, world_radius)
    
    CL_count = csv.sum()
    rate = CL_count / time
    error = np.sqrt(CL_count) / time
    rate = rate.to_numpy()
    ###draw coin, win=20
    rate_win20 = rate#[28:35]
    error_win20 = error#[28:35]

    csv = pd.read_csv(save_dir + "250m.csv")
    event_num = num_per_run * len(csv)
    time = time_estimate(event_num, 250)

    cl1_count = (csv.sum())[0]
    rate_win20[0] = cl1_count / time
    error_win20[0] = np.sqrt(cl1_count) / time

    # pc_cl.ax.step(1, 71863.11901504689, marker='.', color='b', label='world radius = 250m')
    # pc_cl.ax.step(k40coin_num, rate_win20, label='K40', where='mid', color='cornflowerblue')
    pc_cl.ax.plot(k40coin_num, rate_win20, drawstyle='steps-mid', color='cornflowerblue',label='K40')
    print(rate_win20)
    pc_cl.ax.fill_between(k40coin_num, rate_win20,step='mid', color='cornflowerblue', alpha=0.3)

    # # pc_cl.ax.errorbar(1, 71863.11901504689, yerr=12396.756747187115, fmt='.', color='b')
    pc_cl.ax.errorbar(k40coin_num, rate_win20, yerr=error_win20, fmt='.', color='cornflowerblue')

def muon(pc_cl):
    hits_0_16 = pd.read_csv('../domhits0-16.csv',header=None)
    hits_17_57 = pd.read_csv('../domhits17-57.csv',header=None)
    hits_58_123 = pd.read_csv('../domhits58-123.csv',header=None)
    hits = pd.concat([hits_0_16,hits_17_57,hits_58_123])
    # hits = hits_0_16
    dom_num = 124*21
    hits[20] = hits[20] * np.pi * 350**2 / dom_num
    # hits = hits[(hits.index % 21).isin([0,1,2,3,4,5])]

    hits_ = []
    for i in range(20):
        hits_.append((hits[i]*hits[20]**2).sum())
    yerr = np.sqrt(hits_)

    for i in range(20):
            hits[i] = hits[i]*hits[20]
    hits = hits.drop(20,axis=1)
    # yerr = hits.std()
    print(yerr)
    
    hits = hits.sum()

    cl = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    # # pc.ax.plot(cl,hits_2700.sum()*np.pi*350**2/string_num, drawstyle='steps-mid', label='2700m')
    pc_cl.ax.plot(cl,hits, drawstyle='steps-mid',color='orange',label='muon')
    pc_cl.ax.fill_between(cl,hits,step='mid', color='orange', alpha=0.3)
    pc_cl.ax.errorbar(cl,hits,yerr=yerr,fmt='.',color='orange')
    

if __name__ == "__main__":
    
    # pc_cl = PlotContainer(xlabel='coincidence level (CL)', ylabel=r'rate [Hz]', logy=True,
    #                       title="coincidence rate of dom at 2.7km-2.85km depth", figname='all_2.7-2.85km.jpg')

    pc_cl = PlotContainer(xlabel='coincidence level (CL)', ylabel=r'rate [Hz]', logy=True,
                           title="average coincidence rate of doms", figname='all.jpg')


    k40(pc_cl)
    muon(pc_cl)

    # pc_cl.ax.set_xticks(range(1,21,4))
    # pc_cl.ax.set_yticks(np.logspace(-3, 5, 5))
    
    pc_cl.ax.set_xticks(range(1,21,4))
    pc_cl.apply_settings()
    pc_cl.ax.yaxis.set_major_locator(ticker.LogLocator(base=10,numticks=10))
    pc_cl.ax.yaxis.set_minor_locator(ticker.LogLocator(base=10,subs=np.arange(2,10)*.1,numticks=10))
    
    # pc_cl.ax.set_yscale('log')
    # pc_cl.ax.set_yticks(np.logspace(-3, 5, 9))
    

    pc_cl.ax.legend()
    
    pc_cl.savefig()