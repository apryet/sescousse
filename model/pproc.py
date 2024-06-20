import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.ticker import LinearLocator, FormatStrFormatter
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

ax1.bar(swb.index,swb.D,color='darkgrey',label='Recharge')
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
fig.savefig('swb_vars.pdf',dpi=300)



# -------------------------------------
#--- sim obs records  
# -------------------------------------

fig,axs=plt.subplots(2,3,figsize=(10,6),sharex=True)

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

fig.tight_layout()

fig.savefig(os.path.join(sim_dir,'fig','sim_obs_ts.pdf'),dpi=300)


# -------------------------------------
#--- drn records 
# -------------------------------------

fig,axs = plt.subplots(3,1,sharex=True, figsize=(10,6)) #A4 paper size

ax0, ax1, ax2 = axs

ax0.bar(clim.index,clim.D,color='darkgrey',label='Recharge')
ax0.set_ylabel('mm/d')

ax1.bar(hsim.date, drnf_records, color='tan',label='Drainage',alpha=0.8)
ax1.set_ylabel('mm/d')

for loc in ['PS1','PS2','PS3']:
    hsim.loc[:,loc].plot(ax=ax2,label=loc)

ax2.set_ylabel('m NGF')

ax2.set_xlim(hsim.index.min(),hsim.index.max())

fig.align_ylabels()
lgd = [ax.legend(loc='upper left') for ax in [ax0,ax1,ax2]]
fig.tight_layout()
fig.savefig(os.path.join('fig','sim_records.pdf'),dpi=300)


# -------------------------------------
#--- 2D map plot 
# -------------------------------------

# load simulation data 
sim = flopy.mf6.MFSimulation.load(sim_ws='.')
ml = sim.get_model()
nper = sim.tdis.nper.data
hds = ml.output.head().get_alldata()

# 2D map
for i,n in enumerate(range(0,nper,5)):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    modelmap = flopy.plot.PlotMapView(model=ml, ax=ax)
    ax.set_xlim(385600., 387500.)
    ax.set_title('Piézométrie du ' + hsim.date[n].strftime("%d-%m-%Y"))
    pa = modelmap.plot_array(hds[n,:,:],vmin=20,vmax=22)
    cb = plt.colorbar(pa, shrink=0.5)
    cb.set_label('m NGF')
    fig.savefig(os.path.join(sim_dir,'fig',f'h_{i}.png'))

'''
convert 'h_%d.png[0-54]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hmap.gif
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

'''
convert 'hsurf_%d.png[0-54]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hsurf.gif
'''


#

# zone budget 
zarr = np.zeros(ml.modelgrid.shape, dtype=int)
idx = (X > 385600) & (X < 387500)
idx = idx.reshape(zarr.shape)
zarr[idx]=1
zonbud = ml.output.zonebudget(zarr)

cbc = ml.output.budget()

drnb = ml.output.budget().get_data(text='DRN')
drn_surf = 1071422 # m2
drnf_records = -1*np.array([ drnb[i]['q'].sum() for i in range(nper)])/drn_surf*1000*86400 # m3/s to mm/d


# export simulated heads to shapefile 
flopy.export.shapefile_utils.write_grid_shapefile('heads.shp',
                                                  ml.modelgrid,
                                                  {'hmin':hds.min(axis=0)[0,:,:],
                                                   'hmax':hds.max(axis=1)[0,:,:]},
                                                  crs='epsg:2154',
                                                  )



