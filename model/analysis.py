import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy


# load simulation data 
sim = flopy.mf6.MFSimulation.load(sim_ws='.')
ml = sim.get_model()
nper = sim.tdis.nper.data

# set output control and re-reun model 
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}
sim.write_simulation()
#success, buff = sim.run_simulation(report=True)

# load idomain raster
izone_file = os.path.join('..','..','gis','idomain.tif') 
izone_rast = flopy.utils.Raster.load(izone_file).get_array(1)
izone = izone_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')

# load head data 
nper = sim.tdis.nper.data
hdfile = ml.output.head()
hds = hdfile.get_alldata()
times = hdfile.get_times()

# start date
start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())
dates_out  = (start_date + pd.to_timedelta(times,'s')).date

# load top (digital terrain elevation)
top = ml.dis.top.data

# load drain flow 
cbc = ml.output.budget()
drnb = ml.output.budget().get_data(text='DRN')


# -------------------------------------
#---  maps of groundwater head
# -------------------------------------

# 2D map of gw heads 
for n,i in enumerate(range(0,hds.shape[0],5)):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    modelmap = flopy.plot.PlotMapView(model=ml, ax=ax)
    ax.set_xlim(385600., 387500.)
    ax.set_title('Piézométrie du ' + dates_out[i].strftime("%d-%m-%Y"))
    # masked 2D head array
    hds2d = hds[i,0,:,:]
    hds2d[izone < 2]=np.nan
    pa = modelmap.plot_array(hds2d, vmin=20,vmax=21)
    cb = plt.colorbar(pa, shrink=0.5)
    cb.set_label('m NGF')
    fig.savefig(os.path.join('fig',f'h_{n}.png'))
    plt.close()

convert 'h_%d.png[0-54]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hmap.gif


# -------------------------------------
#---  map of groundwater depth
# -------------------------------------

gwd = top - hds

# 2D map of gw depth 
for n,i in enumerate(range(0,hds.shape[0],5)):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    modelmap = flopy.plot.PlotMapView(model=ml, ax=ax)
    ax.set_xlim(385600., 387500.)
    ax.set_title('Profondeur nappe / sol ' + dates_out[i].strftime("%d-%m-%Y"))
    gwd2d = gwd[i,0,:,:]
    gwd2d[izone < 2]=np.nan
    pa = modelmap.plot_array(gwd2d, vmin=0,vmax=2,cmap='viridis')
    cb = plt.colorbar(pa, shrink=0.5)
    cb.set_label('m NGF')
    fig.savefig(os.path.join('fig',f'gwd_{n}.png'))
    plt.close()

convert 'gwd_%d.png[0-54]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither gwdmap.gif

# -------------------------------------
#--- drn records 
# -------------------------------------

drn_surf = 1071422 # m2
drnf_records = -1*np.array([ drnb[i]['q'].sum() for i in range(nper)])/drn_surf*1000*86400 # m3/s to mm/d

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

#
'''
# zone budget 
zarr = np.zeros(ml.modelgrid.shape, dtype=int)
idx = (X > 385600) & (X < 387500)
idx = idx.reshape(zarr.shape)
zarr[idx]=1
zonbud = ml.output.zonebudget(zarr)



# export simulated heads to shapefile 
flopy.export.shapefile_utils.write_grid_shapefile('heads.shp',
                                                  ml.modelgrid,
                                                  {'hmin':hds.min(axis=0)[0,:,:],
                                                   'hmax':hds.max(axis=1)[0,:,:]},
                                                  crs='epsg:2154',
                                                  )


'''
