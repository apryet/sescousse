import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# load simulation data 
sim_dir = '.'
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir,load_only=['dis','tdis'],verbosity_level=0)
start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())
hsim = pd.read_csv(os.path.join(sim_dir,'sescousse.head.csv'),index_col='time')
hsim['date'] = (pd.to_datetime(start_date)+pd.to_timedelta(hsim.index.values.astype(float),'s')).date
hsim.index = hsim.date

# load observations
obs_file = os.path.join(sim_dir,'obs_levels.csv')
hobs = pd.read_csv(obs_file,index_col=0,parse_dates=True)

# load swb data
swb = pd.read_csv(os.path.join('swb_vars.csv'),index_col=0,parse_dates=True)

fig_dir = os.path.join(sim_dir,'fig')

if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

# -----------------------------------------------------------
# soil water budget 
# -----------------------------------------------------------

fig,axs = plt.subplots(4,1,sharex=True, figsize=(10,6)) #A4 paper size

ax0, ax1, ax2, ax3 = axs

# precipitations
ax0.bar(swb.index, swb.P,color='darkblue',label='Precipitations')
ax0.set_ylabel('P [mm/d]')
ax0.legend(loc='upper right')

# potential, soil evapotranspiration and root uptake 
ax1.bar(swb.index, swb.PET, color='grey', label='PET')
ax1.bar(swb.index, swb.RU, color='darkgoldenrod', label='RU',alpha=0.7)
ax1.bar(swb.index, swb['T'], color='forestgreen', label='CSET',alpha=0.7)
ax1.set_ylabel('(P)ET mm/d')
ax1.legend(loc='upper left')

ax2.bar(swb.index,swb.R,color='darkgrey',label='Recharge')
ax2.set_ylabel('mm/d')
lhr,llr = ax2.get_legend_handles_labels()

twax2 = ax2.twinx()
twax2.plot(swb.index,swb.S,color='blue',label='Soil water storage')
twax2.set_ylabel('mm')
lhs,lls = twax2.get_legend_handles_labels()

ax2.legend(lhr+lhs,llr+lls,loc='upper right')

for pid in ['PS1','PS2','PS3']:
    hobs[pid].plot(ax=ax3,label=pid)

ax3.set_ylabel('m NGF')

ax3.set_xlim(swb.index.min(),swb.index.max())

ax3.legend(loc='upper left')

fig.align_ylabels()
fig.tight_layout()
fig.savefig(os.path.join('fig','swb_vars.pdf'),dpi=300)

# full period water budget
swb_sum = swb.iloc[1:].sum(axis=0)
P, PET, T, D, RU, R = [swb_sum[var] for var in ['P','PET','T','D','RU','R']]

DS = swb.S[-1] - swb.S[0]
print(f'Total precipitation :{P:.0f} mm')
print(f'Potential evapotranspiration : {PET:.0f} mm')
print(f'Recharge : {R:.0f} mm')
print(f'Transpiration from soil : {T:.0f} mm')
print(f'Transpiration from aquifer : {RU:.0f}')
print(f'T+RU : {T+RU:.0f} mm')
print(f'P - (T+R) : {P-(T+R):.0f} mm')
print(f'Soil water storage variation:{DS:.0f} mm')

# -------------------------------------
#--- sim obs records  
# -------------------------------------

fig,axs=plt.subplots(2,3,figsize=(10,6),sharex=True,sharey=True)

for ax,loc in zip(axs.ravel(),hsim.columns):
    # plot sim
    hsim[loc].plot(ax=ax, color='darkred',
                   legend=False, label=loc)
    if loc.startswith('FS'):
        hobs_ss = hobs.loc[hobs.index>pd.to_datetime('2023-11-01'),loc]
    else : 
        hobs_ss = hobs[loc]
    # plot obs 
    hobs_ss.plot(ax=ax,marker='+',
                                 legend=False,
                                 ls='',color='darkblue')
    ax.set_xlim(hsim.date.min(),hsim.date.max())
    ax.set_title(loc)
    #ax.set_ylim(19.,22.4)

axs[0,0].axhline(20.86,ls='--',color='grey')
axs[0,1].axhline(20.77,ls='--',color='grey')
axs[0,2].axhline(21.15,ls='--',color='grey')

axs[0,0].set_ylabel('m NGF')
axs[1,0].set_ylabel('m NGF')

fig.tight_layout()
fig.savefig(os.path.join(sim_dir,'fig','sim_obs_ts.pdf'),dpi=300)

