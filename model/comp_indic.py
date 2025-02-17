import os
import numpy as np
import pandas as pd
import matplotlib as mpl
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
    # alert thresholds (depths)
    depth_w = 0.4 # water excess 
    depth_d = 1.2 # water deficit
    # load mf model and simulated heads
    # skip all packages but tdis to save time (otherwise takes ages!)
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir,load_only=['tdis'],verbosity_level=0)
    ml = sim.get_model()
    nper = sim.tdis.nper.data
    headfile = flopy.utils.binaryfile.HeadFile(os.path.join(sim_dir,'sescousse.hds'))
    hds = headfile.get_alldata()
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

# ---------------------------------------------------------
# ----   compute  indicators / calibration period    ------
# ---------------------------------------------------------

# comparison sim dirs 
dirs = {'cal':'master_glm',
        'nodrn':'drn0',
        'drn40':'drn40',
        'drn110':'drn110',
        'drn150':'drn150',
        'drn200':'drn200',
        'drn300':'drn300'
        }


# compute indicators 
indics={}
for l,d in dirs.items():
    indics[l] = get_indic(d)

# -----------------------------------------------------------------------------------
# --    define a list of sim dirs for comparison and a label for this set    --------
# -----------------------------------------------------------------------------------

# historical period  
dirs = [ f'histo_drn{dcm}' for dcm in [0,40,110,150,200,300] ]
label='histo_dd'

# alternative drainage networks
dirs = ['histo_drn40','histo_drn110',\
        'histo_drn_alt_ld_40','histo_drn_alt_ld_110',\
        'histo_drn_alt_hd_40','histo_drn_alt_hd_110']

label='histo_alt'

# prospective periods
dirs = [f'prosp_{i:02d}' for i in range(1,12)]
label='prosp'

# ------------------------------------------------------------------
# --    compute indicators    --------
# ------------------------------------------------------------------

from multiprocessing import Manager, Pool

def store_indic(args):
    d, indics, lock = args
    print(f'Starting store_dic for {d}')
    df=get_indic(d)
    with lock: # ensure only one process modifies the dictionary at a time
        print(f'Storing indics for {d}')
        indics[d] = df 

manager = Manager()
indics= manager.dict()
lock = manager.Lock()  

with Pool(processes=2) as pool:  # Limit the number of parallel processes to 2
    pool.map(store_indic, [(d, indics, lock) for d in dirs])


# ---------------------------------------------------------
# ----   compute  drainage rate     --------------
# ---------------------------------------------------------

zone_surf = 1071422 # m2
drnobs = {}

for d in dirs:
    if ('drn0' in d) : continue # skip simulations without drain package  
    df = pd.read_csv(os.path.join(d,'sescousse.drn.obs.output.csv'),index_col=0)
    df.index = pd.to_datetime(indics[d].index)
    df['flow']= df.sum(axis=1)*(-1/zone_surf*1000*86400) # m3/s to mm/d
    drnobs[d]=df


# ---------------------------------------------------------
# ----    plots    --------------
# ---------------------------------------------------------

'''
colors = {'cal':'k','nodrn':'darkblue',
          'drn40':'darkgreen','drn110':'tan','drn150':'orange','drn200':'darkred','drn300':'purple'}
'''

cmap = mpl.colormaps['tab10']
colors =cmap(np.linspace(0, 1, len(dirs)))
colors = ['darkblue','darkgreen','tan','orange','darkred','purple']
colors_dic = {d:c for d,c in zip(dirs,colors)}


# --- daily records of indicators 
fig,ax =plt.subplots(1,1,figsize=(dbcol_width,0.5*dbcol_width))

for d in dirs:
    df = indics[d]
    dates_out, w_records, d_records = df.index,df.w,df.d
    w_records = w_records.round(2)
    d_records = d_records.round(2)
    w_records[w_records==0]=np.nan
    d_records[d_records==0]=np.nan
    ax.plot(dates_out,w_records,color=colors_dic[d],ls='-',label=d,alpha=0.8)
    ax.plot(dates_out,d_records,color=colors_dic[d],ls='-',alpha=0.8)

ax.set_ylabel('Deficit / Excess [m]')
ax.legend()
ax.axhline(0,linewidth=1,color='k',linestyle='--')
fig.tight_layout()
fig.savefig(os.path.join('fig',f'{label}_indics_records.pdf'),dpi=300)


# --- cumulated annual indicator values 

nyears = round((indics[dirs[0]].index[-1]-indics[dirs[1]].index[0]).days/365.2425)
indics_cum = pd.DataFrame({
    'wcum':[ indics[d]['w'].sum() for d in dirs],
    'dcum':[ indics[d]['d'].sum() for d in dirs]
    }, index = dirs).div(nyears)

fig,ax = plt.subplots(1,1,figsize=(mdcol_width,mdcol_width))
indics_cum.plot(ax=ax,kind='bar',color=['darkblue','tan'])
ax.set_ylabel('Deficit / Excess [m$\\times$j/an]')
fig.tight_layout()
fig.savefig(os.path.join('fig',f'{label}_indics_cum.pdf'),dpi=300)

# --- daily records of drainage rate 
fig,ax =plt.subplots(1,1,figsize=(dbcol_width,0.5*dbcol_width))

for d in drnobs.keys():
    df = drnobs[d]
    ax.plot(df.index, df.flow, color=colors_dic[d],label=d,alpha=0.8)
    ax.legend()

ax.set_ylabel('Drainage rate [mm/d]')
fig.savefig(os.path.join('fig',f'{label}_drn_records.pdf'),dpi=300)

# save into single df 
drnobs_df = pd.concat(drnobs.values(), axis=1, keys=drnobs.keys())
dd = drnobs_df.loc[:,(slice(None),'flow')]
dd.columns = dd.columns.get_level_values(0)
dd.to_excel(os.path.join('fig',f'{label}_sim_histo_drn.xlsx'))

# --- cumulated annual drainage values 

drnobs_cum = pd.DataFrame({d:drnobs[d].flow.groupby(pd.Grouper(freq='12MS',origin=drnobs[d].index[0])).sum() for d in drnobs.keys()})

fig,ax = plt.subplots(1,1,figsize=(mdcol_width,0.5*mdcol_width))
drnobs_cum.plot(ax=ax,color=[colors_dic[d] for d in drnobs_cum.columns])
ax.set_ylabel('Drainage rate [mm/year]')
ax.set_xticklabels(list(drnobs_cum.index.year))
fig.tight_layout()
fig.savefig(os.path.join('fig',f'{label}_drn_cum.pdf'),dpi=300)

