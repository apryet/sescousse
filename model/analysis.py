import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.dates as mdates
import flopy

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54@
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# load simulation data 
sim = flopy.mf6.MFSimulation.load(sim_ws='.')
ml = sim.get_model()
ml.modelgrid.set_coord_info(crs=2154)
nper = sim.tdis.nper.data
delr = ml.dis.delr[0]
delc = ml.dis.delc[0]

'''
# set output control and re-reun model 
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}
sim.write_simulation()
success, buff = sim.run_simulation(report=True)
'''
# load idomain raster
idomain_file = os.path.join('..','..','gis','idomain.tif') 
idomain_rast = flopy.utils.Raster.load(idomain_file)
idomain = idomain_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
idomain_3d = np.stack([idomain]*nper) # transient time domain 

# load dtm raster
dtm_file = os.path.join('..','..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)
dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')

# load top (digital terrain elevation)
top = ml.dis.top.data

# load swb 
swb = pd.read_csv(os.path.join('swb_vars.csv'),index_col=0,parse_dates=True)

# load cbc 
cbc = ml.output.budget()

# load head data 
nper = sim.tdis.nper.data
hdfile = ml.output.head()
hdfile.mg = ml.modelgrid
hds = hdfile.get_alldata()
times = hdfile.get_times()

# start date
start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())
end_date = pd.to_datetime(start_date+ pd.to_timedelta(nper-1,'d'))
dates_out  = pd.date_range(start_date,end_date).date

# load simulated heads
hsim = pd.read_csv(os.path.join('sescousse.head.csv'),index_col='time')
hsim['date'] = (pd.to_datetime(start_date)+pd.to_timedelta(hsim.index.values.astype(float),'s')).date
hsim.index = hsim.date

# load observations
obs_file = os.path.join('obs_levels.csv')
hobs = pd.read_csv(obs_file,index_col=0,parse_dates=True)

'''
# -------------------------------------
#--- drn records 
# -------------------------------------

# load drain flow 
#zone_surf = delr*delc*(idomain>1).sum()
zone_surf = 1071422
drnobs = pd.read_csv('sescousse.drn.obs.output.csv',index_col=0)
drnobs.index = dates_out
drnobs['flow']= drnobs.sum(axis=1)*(-1/zone_surf*1000*86400) # m3/s to mm/d

drnb = ml.output.budget().get_data(text='DRN')
drnf_records = -1*np.array([ drnb[i]['q'].sum() for i in range(nper)])/zone_surf*1000*86400 # m3/s to mm/d

fig,axs = plt.subplots(3,1,sharex=True, figsize=(10,6)) #A4 paper size
ax0, ax1, ax2 = axs
ax0.bar(swb.index,swb.R,color='darkgrey',label='Recharge')
ax0.set_ylabel('mm/d')
#ax1.bar(dates_out, drnf_records, color='orange',label='Drainage CBC',alpha=0.8)
ax1.bar(dates_out, drnobs.flow, color='tan',label='Drainage',alpha=0.8)
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
#--- evt records 
# -------------------------------------

# load actual evt from cbc 
evtb = ml.output.budget().get_data(text='EVT',full3D=True) # so slow !
evtba = np.stack(evtb)[:,0,:,:]

# mask out of area of interest 
mevtba = np.ma.masked_where(idomain_3d<3, evtba)

# time series of mean evt
evt_records = -1*np.array(mevtba.mean(axis=(1,2))/(delr*delc)*86400*1000) # m3/s to mm/d

swb['ARU'] = evt_records

# potential evapotranspiration at weekly 

# aggregate by 7-days 
wswb = swb.groupby(pd.Grouper(freq='7D')).sum()
bwdth = 5

fig,axs = plt.subplots(3,1,sharex=True, figsize=(dbcol_width,dbcol_width)) 
ax0, ax1, ax2 = axs

# precipitations
ax0.bar(wswb.index, wswb.P, width=bwdth,color='darkblue',label='Précipitations')
ax0.set_ylabel('P [mm/semaine]')
ax0.legend(loc='upper right')

# potential evapotranspiration
ax1.bar(wswb.index, wswb.PET, width=bwdth, color='tan', label='Evapotranspiration potentielle')
ax1.bar(wswb.index, wswb['T'], width=bwdth, color='olive', label='Evapotranspiration depuis le sol')
ax1.bar(wswb.index, wswb['ARU'], width=bwdth, bottom = wswb['T'], color='darkgreen', label='Prélèvement racinaire en aquifère')
ax1.set_ylabel('(P)ET [mm/semaine]')
ax1.legend(loc='upper left')

# recharge
ax2.bar(wswb.index, wswb.R, width=bwdth, color='darkgrey',label='Recharge de l\'aquifère')
ax2.set_ylabel('Recharge [mm/semaine]')
ax2.legend()

fig.tight_layout()
fig.savefig(os.path.join('fig','wswb_pet_et_ru.pdf'),dpi=300)
'''
# -------------------------------------
#--- indicators
# -------------------------------------
 
# critical depths
depth_w = 0.4 # water excess 
depth_d = 1.5 # water stress

# critical levels 
z_w = dtm - depth_w
z_d = dtm - depth_d

# --- plot head records at PS1 with critical depths 

fig,ax =plt.subplots(1,1,figsize=(dbcol_width,dbcol_width/2))

# get row, col of observation point PS1 
obs_df = pd.DataFrame(ml.obs[0].continuous.data['sescousse.head.csv'])
obs_df = obs_df.set_index('obsname')
i = obs_df.loc['ps1','id'][1] - 1 # 0-based row
j = obs_df.loc['ps1','id'][2] - 1 # 0-based col

ps1_terrain_level = dtm[i,j]
ps1_w_level = z_w[i,j]
ps1_d_level = z_d[i,j]

# gw level
ax.plot(dates_out,hsim.PS1,label='Groundwater level at PS1',color='k')
# surface level
ax.axhline(ps1_terrain_level,color='darkgrey',ls='--')
ax.text(dates_out[-50],ps1_terrain_level+0.05,'Surface level',color='darkgrey')
# wet critical level
ax.axhline(ps1_w_level,color='darkblue',ls=':')
ax.text(dates_out[-50],ps1_w_level+0.05,'Wet depth',color='darkblue')
# dry critical level 
ax.axhline(ps1_d_level,color='darkorange',ls=':')
ax.text(dates_out[-50],ps1_d_level+0.05,'Dry depth',color='darkorange')

ax.legend(loc='upper left')
ymin,ymax=ax.get_ylim()
ax.set_ylim(ymin,ymax+0.20)
ax.set_ylabel('Elevation [m NGF]')

fig.savefig(os.path.join('fig','critical_levels_PS1.pdf'),dpi=300)

'''
# --- plot indicator records 

# masked head array over area of interest (idomain=3)
mhds = np.ma.masked_where(idomain_3d<3, hds[:,0,:,:])

# spatially averaged dtm 
mdtm = np.ma.masked_where(idomain<2, dtm)
zm = mdtm.mean()

# spatially averaged gw levels 
hm = mhds.mean(axis=(1,2))

# records of spatially averaged water excess/stress
ncells = (idomain==3).sum()
w_records = ((mhds-z_w)*(mhds>z_w)).sum(axis=(1,2))/ncells
d_records = ((mhds-z_d)*(mhds<z_d)).sum(axis=(1,2))/ncells


# cumulated water excess/stress 
w = w_records.sum()
d = d_records.sum()

fig,ax =plt.subplots(1,1,figsize=(dbcol_width,dbcol_width/2))

ax.plot(dates_out,w_records,color='darkblue',label='Water excess')
ax.plot(dates_out,d_records,color='darkorange',label='Water deficit')
ax.set_ylabel('Water excess / deficit [m]')
ax.legend()
fig.tight_layout()

fig.savefig(os.path.join('fig','indicators_records.pdf'),dpi=300)


# -------------------------------------
#---  save heads to shapefiles
# -------------------------------------

gis_dir = os.path.join('gis')

if not os.path.exists(gis_dir):
    os.mkdir(gis_dir)

mhds = hds.mean(axis=(1,2,3))
iper_dry = mhds.argmin()
iper_wet = mhds.argmax()

head_shpfile = os.path.join(gis_dir,'heads_ss.shp')
hdfile.to_shapefile(head_shpfile,kstpkper=(0,0))

head_shpfile = os.path.join(gis_dir,'heads_dry.shp')
hdfile.to_shapefile(head_shpfile,kstpkper=(0,iper_dry))

head_shpfile = os.path.join(gis_dir,'heads_wet.shp')
hdfile.to_shapefile(head_shpfile,kstpkper=(0,iper_wet))

# ------------------------------------
#---  maps
# -------------------------------------

# -------------------------------------
#---  maps of groundwater head
# -------------------------------------

hmaps_dir = os.path.join('fig','hmaps')

if not os.path.exists(hmaps_dir):
    os.mkdir(hmaps_dir)

# 2D map of gw heads 
for n,i in enumerate(range(0,hds.shape[0],1)):
    fig,axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]}, figsize=(7,7)) 
    # records with time line 
    ax0, ax1 = axs
    ax0.bar(swb.index,swb.P,color='darkblue',label='Precipitations')
    ax0.set_ylabel('mm/d')
    ax0.axvline(dates_out[i],color='k',linewidth=2,linestyle='--')
    twax0 = ax0.twinx()
    twax0.plot(hobs.index,hobs.FS4,color='darkred',label='FS4')
    twax0.set_ylabel('m NGF')
    twax0.set_xlim(swb.index.min(),swb.index.max())
    twax0.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax0.xaxis.get_major_locator()))
    # map
    modelmap = flopy.plot.PlotMapView(model=ml, ax=ax1)
    ax1.set_xlim(385600., 387500.)
    #ax1.set_title('Piézométrie du ' + dates_out[i].strftime("%d-%m-%Y"))
    # masked 2D head array
    hds2d = hds[i,0,:,:]
    hds2d[idomain < 2]=np.nan
    pa = modelmap.plot_array(hds2d, vmin=20.5,vmax=22)
    cb = plt.colorbar(pa, shrink=0.5)
    cb.set_label('m NGF')
    ax1.set_aspect(1)
    #ax1.set_xlabels([])
    #ax1.set_ylabels([])
    fig.tight_layout()
    fig.savefig(os.path.join(hmaps_dir,f'h_{n}.png'))
    plt.close()

#convert 'h_%d.png[0-50]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hmap.gif

# -------------------------------------
#---  map of groundwater depth
# -------------------------------------

gwdmaps_dir = os.path.join('fig','gwdmaps')

if not os.path.exists(gwdmaps_dir):
    os.mkdir(gwdmaps_dir)

gwd = dtm - hds

# 2D map of gw depth 
for n,i in enumerate(range(0,hds.shape[0],1)):
    fig,axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]}, figsize=(9,9)) 
    # records with time line 
    ax0, ax1 = axs
    ax0.bar(swb.index,swb.P,color='darkblue',label='Precipitations')
    ax0.set_ylabel('mm/d')
    ax0.axvline(dates_out[i],color='k',linewidth=2,linestyle='--')
    twax0 = ax0.twinx()
    twax0.plot(hobs.index,hobs.FS4,color='darkred',label='FS4')
    twax0.set_ylabel('m NGF')
    twax0.set_xlim(swb.index.min(),swb.index.max())
    twax0.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax0.xaxis.get_major_locator()))
    # map
    modelmap = flopy.plot.PlotMapView(model=ml, ax=ax1)
    ax1.set_xlim(385600., 387500.)
    ax1.set_title('Profondeur nappe / sol ' + dates_out[i].strftime("%d-%m-%Y"))
    gwd2d = gwd[i,0,:,:]
    gwd2d[idomain < 2]=np.nan
    pa = modelmap.plot_array(gwd2d, vmin=0,vmax=2,cmap='viridis')
    cb = plt.colorbar(pa, shrink=0.5)
    cb.set_label('Groundwater depth')
    ax1.set_aspect(1)    
    fig.savefig(os.path.join(gwdmaps_dir,f'gwd_{n}.png'))
    plt.close()

#convert 'gwd_%d.png[0-50]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither gwdmap.gif
'''

# -------------------------------------
#--- 3D surface plot 
# -------------------------------------

'''
X = ml.modelgrid.xcellcenters
Y = ml.modelgrid.ycellcenters

surfmaps_dir = os.path.join('fig','surfmaps')

if not os.path.exists(surfmaps_dir):
    os.mkdir(surfmaps_dir)

fig, ax = plt.subplots(1,1,figsize=(12,12),subplot_kw={"projection": "3d"})
ax.view_init(elev=15., azim=-148)
for i,n in enumerate(range(0,nper,7)):
    Z = hds[n,:,:][0,:,:]
    if i>0 : surf.remove()
    surf = ax.plot_surface(X, Y, Z, 
                    rstride=1, cstride=1, 
                    cmap=cm.viridis,
                           vmin=20.5,vmax=22,
                    linewidth=0, antialiased=False)
    ax.set_xlim(385800., 387500.)
    ax.set_zlim(19,23)
    ax.set_title('Piézométrie du ' + hsim.date[n].strftime("%d-%m-%Y"))
    if i==0 : 
        cbar = fig.colorbar(surf, shrink=0.3, aspect=10)
        cbar.set_label('h [m NGF]')
    fig.savefig(os.path.join(surfmaps_dir,f'hsurf_{i}.png'),dpi=128)

#convert 'hsurf_%d.png[0-50]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hsurf.gif

# -------------------------------------
#---  x-section through PS1
# -------------------------------------

xsects_dir = os.path.join('fig','xsects')

if not os.path.exists(xsects_dir):
    os.mkdir(xsects_dir)


# get row, col of observation point PS1 
obs_df = pd.DataFrame(ml.obs[0].continuous.data['sescousse.head.csv'])
obs_df = obs_df.set_index('obsname')
ps1_row = obs_df.loc['ps1','id'][1] - 1 # 0-based col

hmin,hmax=19.5,22

for n,i in enumerate(range(0,hds.shape[0],1)):
    fig,axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]}, figsize=(7,7)) 
    # records with time line 
    ax0, ax1 = axs
    ax0.bar(swb.index,swb.P,color='darkblue',label='Precipitations')
    ax0.set_ylabel('mm/d')
    ax0.axvline(dates_out[i],color='k',linewidth=2,linestyle='--')
    twax0 = ax0.twinx()
    twax0.plot(hobs.index,hobs.FS4,color='darkred',label='FS4')
    twax0.set_ylabel('m NGF')
    twax0.set_xlim(swb.index.min(),swb.index.max())
    twax0.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    # map
    hds2d = hds[i,0,:,:]
    xsect = flopy.plot.PlotCrossSection(model=ml,ax=ax1,line={"Row": ps1_row})
    pc = xsect.plot_array(hds2d, head=hds2d, alpha=1)
    bc=  xsect.plot_bc('drn',color='tan',alpha=0.7)
    bc=  xsect.plot_bc('riv',color='blue',alpha=0.7)
    ax1.set_xlim(190,1040)
    ax1.set_ylim(hmin,hmax)
    linecollection = xsect.plot_grid(alpha=0.5)
    cb = plt.colorbar(pc, shrink=0.75)
    cb.set_label('m NGF')
    fig.savefig(os.path.join(xsects_dir,f'xsect_{n}.png'))
    plt.close()

#convert 'xsect_%d.png[0-352]' -scale 1066x800 -delay 20 -coalesce -layers Optimize -fuzz 2% +dither hmap.gif


'''
