import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# load idomain raster
idomain_file = os.path.join('..','gis','idomain.tif') 
idomain_rast = flopy.utils.Raster.load(idomain_file)

# load dtm raster
dtm_file = os.path.join('..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)

# -------------------------------------
# ----  indicator definitions     -----
# -------------------------------------

def get_indic(sim_dir):
    # critical depths
    depth_w = 0.4 # water excess 
    depth_d = 1.5 # water deficit
    # load mf model and simulated heads
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir)
    ml = sim.get_model()
    nper = sim.tdis.nper.data 
    hds = ml.output.head().get_alldata()
    # resample idomain and dtm 
    idomain = idomain_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
    idomain_3d = np.stack([idomain]*nper) # transient idomain 
    dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
    # number of cells in domain of interest (idomain==3)
    ncells = (idomain==3).sum()
    # critical levels 
    z_w = dtm - depth_w # m NGF
    z_d = dtm - depth_d # m NGF
    # dates out 
    start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())
    end_date = pd.to_datetime(start_date+ pd.to_timedelta(nper-1,'d'))
    dates_out  = pd.date_range(start_date,end_date).date
    # masked head array out of area of interest (idomain==3)
    mhds = np.ma.masked_where(idomain_3d<3, hds[:,0,:,:])
    # records of spatially averaged water excess/stress
    w_records = ((mhds-z_w)*(mhds>z_w)).sum(axis=(1,2))/ncells*1000
    d_records = ((mhds-z_d)*(mhds<z_d)).sum(axis=(1,2))/ncells*1000
    # aggregate in df
    df = pd.DataFrame({'w':w_records,'d':d_records},index=dates_out)
    return(df)

# -------------------------------------
# ----    indicators     --------------
# -------------------------------------

# comparison sim dirs 
dirs = {'cal':'master_glm',
        'nodrn':'nodrn',
        'drn110':'drn110',
        'drn40':'drn40'}

# compute indicators 
indics={}
for l,d in dirs.items():
    indics[l] = get_indic(d)


# plot 
lss = {'cal':'-','nodrn':':','drn110':'--','drn40':'-.'}

fig,ax =plt.subplots(1,1,figsize=(dbcol_width,mdcol_width))

dirs = {'drn110':'drn110',
        'drn40':'drn40'}


for l,d in dirs.items():
    df = indics[l]
    dates_out, w_records, d_records = df.index,df.w,df.d
    ax.plot(dates_out,w_records,color='darkblue',ls=lss[l],label=f'{l} Water excess',alpha=0.8)
    ax.plot(dates_out,d_records,color='darkorange',ls=lss[l],label=f'{l} Water deficit',alpha=0.8)
    ax.set_ylabel('Water excess / deficit [mm]')
    ax.legend()

#ax.set_ylim(-5,5)

fig.savefig(os.path.join('fig','indics_records.pdf'),dpi=300)


indics_cum = pd.DataFrame({'wcum':[ df['w'].sum() for df in indics.values()],
'dcum':[ df['d'].sum() for df in indics.values()]},
             index = indics.keys())


fig,ax = plt.subplots(1,1,figsize=(mdcol_width,mdcol_width))
indics_cum.plot(ax=ax,kind='bar')
ax.set_ylabel('Hauteur d\'eau cumulée [mm]')
fig.tight_layout()
fig.savefig(os.path.join('fig','indics_cum.pdf'),dpi=300)


