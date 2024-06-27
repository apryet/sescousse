import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

# load simulation data 
sim_dir = '.'
start_date = pd.to_datetime('2023-07-12')
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

fig,axs = plt.subplots(3,1,sharex=True, figsize=(10,6)) #A4 paper size

ax0, ax1, ax2 = axs

ax0.bar(swb.index, swb.PET, color='tan', label='PET')
ax0.bar(swb.index, swb.RU, color='darkgreen', label='Groundwater uptake')
ax0.set_ylabel('(P)ET mm/d')

ax0.legend(loc='upper left')

twax0 = ax0.twinx()
twax0.bar(swb.index, swb.P,color='darkblue',label='Rainfall')
twax0.invert_yaxis()
twax0.set_ylabel('P [mm/d]')
twax0.legend(loc='upper right')

ax1.bar(swb.index,swb.R,color='darkgrey',label='Recharge')
ax1.set_ylabel('mm/d')

ax1.legend(loc='upper left')

twax1 = ax1.twinx()
twax1.plot(swb.index,swb.S,color='blue',label='Soil water storage')
twax1.set_ylabel('mm')

twax1.legend(loc='upper right')

for pid in ['PS1','PS2','PS3']:
    hobs[pid].plot(ax=ax2,label=pid)

ax2.set_ylabel('mm')

ax2.set_xlim(swb.index.min(),swb.index.max())

ax2.legend(loc='upper left')

fig.align_ylabels()
fig.tight_layout()
fig.savefig(os.path.join('fig','swb_vars.pdf'),dpi=300)

# full period water budget
swb_sum = swb.iloc[1:].sum(axis=0)
P, PET, T, D, RU, R = [swb_sum[var] for var in ['P','PET','T','D','RU','R']]

DS = swb.S[-1] - swb.S[0]
print(f'Total precipitation :{P:.0f} mm')
print(f'Potential evapotranspiration : {PET:.0f} mm')
print(f'Recharge : {D:.0f} mm')
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
    ax.set_ylim(19.6,22.4)


axs[0,0].set_ylabel('m NGF')
axs[1,0].set_ylabel('m NGF')
#axs[0,0].axhline(21.46,ls='--',color='grey')
#axs[0,1].axhline(21.37,ls='--',color='grey')
#axs[0,2].axhline(21.75,ls='--',color='grey')

fig.tight_layout()
fig.savefig(os.path.join(sim_dir,'fig','sim_obs_ts.pdf'),dpi=300)


# -------------------------------------
#--- 2D map plot 
# -------------------------------------

# load simulation data 
sim = flopy.mf6.MFSimulation.load(sim_ws='.')

ml = sim.get_model()

'''
oc = ml.get_package('oc')
oc.saverecord = [("HEAD", "ALL")]
oc.head_filerecord = f'{ml.name}.hds'
sim.write_simulation()
'''

nper = sim.tdis.nper.data
hds = ml.output.head().get_alldata()
kstepkper = ml.output.head().get_kstpkper()



'''
# -------------------------------------
#--- 3D surface plot 
# -------------------------------------

X = ml.modelgrid.xcellcenters
Y = ml.modelgrid.ycellcenters

fig, ax = plt.subplots(1,1,figsize=(12,12),subplot_kw={"projection": "3d"})
for i,n in enumerate(range(0,nper,5)):
    Z = hds[n,:,:][0,:,:]
    if i>0 : surf.remove()
    surf = ax.plot_surface(X, Y, Z, 
                    rstride=1, cstride=1, 
                    cmap=cm.viridis,
                    linewidth=0, antialiased=False)
    ax.set_xlim(385600., 387500.)
    ax.set_zlim(18,24)
    ax.set_title('Piézométrie du ' + hsim.date[n].strftime("%d-%m-%Y"))
    if i==0 : 
        cbar = fig.colorbar(surf, shrink=0.3, aspect=10)
        cbar.set_label('h [m NGF]')
    fig.savefig(os.path.join(sim_dir,'fig',f'hsurf_{i}.png'),dpi=128)

#convert 'hsurf_%d.png[0-54]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hsurf.gif



