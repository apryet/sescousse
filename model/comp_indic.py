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
    w_records = ((mhds-z_w)*(mhds>z_w)).sum(axis=(1,2))/ncells
    d_records = ((mhds-z_d)*(mhds<z_d)).sum(axis=(1,2))/ncells
    # aggregate in df
    df = pd.DataFrame({'w':w_records,'d':d_records},index=dates_out)
    return(df)

# -------------------------------------
# ----    indicators     --------------
# -------------------------------------

# comparison sim dirs 
dirs = {'cal':'master_glm',
        'nodrn':'nodrn',
        'drn40':'drn40',
        'drn110':'drn110',
        'drn210':'drn210',
        'drn300':'drn300'
        }


# comparison sim dirs 
dirs = {'nodrn':'nodrn_histo',
        'drn40':'drn40_histo',
        'drn110':'drn110_histo',
        }


# comparison sim dirs 
dirs = {
        'drn110':'drn110_histo'
        }
        


# compute indicators 
indics={}
for l,d in dirs.items():
    indics[l] = get_indic(d)


'''
from multiprocessing import Process

processes = []

for m in range(1,16):
   n = m + 1
   p = Process(target=some_function, args=(m, n))
   p.start()
   processes.append(p)

'''


# plot 
colors = {'cal':'k','nodrn':'darkblue','drn40':'darkgreen','drn110':'red','drn210':'darkred','drn300':'purple'}
lss = {'cal':'-','nodrn':'-','drn40':'-','drn110':'-','drn210':'-','drn300':'-'}


fig,ax =plt.subplots(1,1,figsize=(dbcol_width,0.5*dbcol_width))

for l,d in dirs.items():
    df = indics[l]
    dates_out, w_records, d_records = df.index,df.w,df.d
    ax.plot(dates_out,w_records,color=colors[l],ls=lss[l],label=f'{l}',alpha=0.8)
    ax.plot(dates_out,d_records,color=colors[l],ls=lss[l],alpha=0.8)
    ax.set_ylabel('Deficit / Excess  [m]')
    ax.legend()

#ax.set_ylim(-5,5)
fig.savefig(os.path.join('fig','indics_records.pdf'),dpi=300)


nyears = 1
indics_cum = pd.DataFrame({
    'wcum':[ indics[l]['w'].sum() for l in dirs.keys()],
    'dcum':[ indics[l]['d'].sum() for l in dirs.keys()]
    }, index = list(dirs.keys())).div(nyears)

fig,ax = plt.subplots(1,1,figsize=(mdcol_width,mdcol_width))
indics_cum.plot(ax=ax,kind='bar',color=['darkblue','tan'])
ax.set_ylabel('Deficit / Excess [m$\\times$j/an]')
fig.tight_layout()
fig.savefig(os.path.join('fig','indics_cum.pdf'),dpi=300)


