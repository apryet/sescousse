import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import swb

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# --------------------------------------
# ---- settings     ----------
# --------------------------------------

# steady model dir
steady_dir = 'ml'

# transient model dir
transient_dir = 'ml_transient'

# load steady state model 
sim = flopy.mf6.MFSimulation.load(sim_ws=steady_dir)
ml = sim.get_model()

# set sim start end dates 
start_date = pd.to_datetime('2023-10-15').date()
end_date = pd.to_datetime('2024-10-01').date()

sim_dates = pd.date_range(start_date,end_date).date

# --------------------------------------
# --- load and pre-proc input data  
# --------------------------------------

# load dtm raster
dtm_file = os.path.join('..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)
dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')

# --- save observed levels over simulation period
levels_file = os.path.join('..','data','levels_daily.csv')
levels = pd.read_csv(levels_file,header=[0,1],index_col=0,parse_dates=True)
obs_locs = ['FS1','FS2','FS3','FS4','PS1','PS2','PS3']
obs_levels = levels.xs('h',axis=1)[obs_locs]
obs_levels.index = obs_levels.index.date
obs_levels.to_csv(os.path.join('ml_transient','obs_levels.csv'))

# load initial parameter values 
pdata = pd.read_excel(os.path.join('..','data','par.xlsx'),index_col='name')
parvals = pdata['val']

# Weather data 
#stjean_file=os.path.join('..','data','stjean_daily.csv')
#stjean_df = pd.read_csv(stjean_file, header=[0,1], index_col=0, parse_dates=True)

# ET0 from Meteo-France (Mérignac Station)
et0_file=os.path.join('..','data','ET0_daily.csv')
et0_df = pd.read_csv(et0_file,header=[0],sep=',')
et0_df.index = pd.DatetimeIndex(pd.to_datetime(et0_df.date,format='%d/%m/%Y')).date
et0 = et0_df['ET0']

# P from Meteo-France (Mérignac Station)
p_file=os.path.join('..','data','P_daily.csv')
p_df = pd.read_csv(p_file, header=[0],sep=',')
p_df.index = pd.DatetimeIndex(pd.to_datetime(p_df.date,format='%d/%m/%Y')).date
p = p_df['P']

# save to single clim file over simulation period
clim_df = pd.DataFrame({'P':p.loc[sim_dates].values,'ET0':et0.loc[sim_dates].values},index=sim_dates)
clim_df.to_csv(os.path.join(transient_dir,'clim.csv'))

# Piezometric records from ADES
# https://ades.eaufrance.fr/Fiche/PtEau?code=07545X0029/F
ades_file=os.path.join('..','data','BSS001VYWT.csv')
ades_df = pd.read_csv(ades_file,header=[0],sep=',')
ades_df.index = pd.DatetimeIndex(pd.to_datetime(ades_df.date,format='%d/%m/%Y')).date
ades = ades_df['F']
ades = ades.reindex(sim_dates)
ades = ades.interpolate()

# --------------------------------------
# --- setup recharge model 
# --------------------------------------

swb.run_swb(theta_sat = parvals.loc['tsat'],
            D_max= parvals.loc['dmax'],
            Kc = parvals.loc['kc'],
            cwd=transient_dir)

swb = pd.read_csv(os.path.join(transient_dir,'swb_vars.csv'),index_col=0,parse_dates=True)
rech = swb.loc[:,'R']*0.001/86400 # mm/d to m/s # recharge
evt = swb.loc[:,'RU']*0.001/86400 # mm/d to m/s # transpiration (root water uptake)

# no recharge nor transpiration for initial steady state
rech[0] = 0
evt[0] = 0

# --- tdis package 

# set start date and stress periods
tdis = sim.get_package('tdis')
tdis.start_date_time.set_data(start_date.isoformat()) #ISO format !
nper = (end_date - start_date).days+1
tdis.nper = nper
perlen = 86400 # s
tdis.perioddata = [ [perlen,1,1] for _ in range(nper)]

# --- recharge package
rcha = ml.get_package('rcha_0')
rcha.recharge = {i:rech[i] for i in range(nper)}

# --- evt package (root aquifer water uptake)

evt_surf = dtm - 1.20 # evt surface elevation
evt_extdp = 0.60 # extinction depth from evt_surf 

evta = flopy.mf6.ModflowGwfevta(
    ml,
    pname="evt",
    save_flows=True,
    surface = evt_surf, # non-limited evt elevation 
    depth = evt_extdp, # exctinction elevation 
    rate = {i:evt[i] for i in range(nper)}
)

# -- storage package 
sto = flopy.mf6.ModflowGwfsto(
    ml,
    pname="sto",
    save_flows=True,
    iconvert=1,
    ss=1e-6,
    sy=parvals.loc['sy'],
    steady_state={0: True},
    transient={i: True for i in range(1,nper)},
)

# --- drn package 
drn = ml.get_package('drn_0')
drn_rec0 = drn.stress_period_data.get_data()[0]
drn_rec0['cond'] = parvals.loc['cdrn']

drn_records = obs_levels.loc[sim_dates,['FS1','FS2','FS3']]

# gen stress period data 
drn_spd = {}
for i in range(nper):
    drn_rec = drn_rec0.copy()
    # assign elevation based on boundname
    drn_ids = np.char.upper(drn_rec0['boundname'].astype('str'))
    drn_rec['elev'] = drn_records.iloc[i].loc[drn_ids]
    drn_spd[i]=drn_rec

ml.drn.stress_period_data.set_data(drn_spd)

# --- river network 

# river level records at FS4
riv_ref_records = obs_levels.loc[sim_dates,'FS4'].values

# original recarray ModflowGwfriv
riv = ml.get_package('riv_0')
riv_rec0 = riv.stress_period_data.get_data()[0]

# approach for river network simulation 
rivtype = 'riv'
#rivtype = 'drn'

if rivtype == 'drn':
    hfield = 'elev'
    # new ref recarray for ModflowGwfdrn
    riv_rec0 = flopy.mf6.ModflowGwfdrn.stress_period_data.empty(ml,maxbound=riv_rec0.shape[0])[0]
    riv_rec0['cellid'] = riv_rec0['cellid']
    rivpckg = flopy.mf6.ModflowGwfdrn
else :
    hfield = 'stage'
    rivpckg = flopy.mf6.ModflowGwfriv

# set conductance value from par file 
riv_rec0['cond'] = parvals.loc['criv']

# reference value at FS4 for steady state simulation
riv_ref_value = 20.66 # hriv(fs4,ss)

# gen stress period data 
riv_spd = {}
for i in range(nper):
    riv_rec = riv_rec0.copy()
    # h_riv = hriv(ss) + (hriv(fs4,t)-hriv(fs4,ss)
    riv_rec[hfield] = riv_rec0[hfield]+ (riv_ref_records[i]-riv_ref_value)
    riv_spd[i]=riv_rec

# remove former package 
ml.remove_package('riv')
# create new package 
riv = rivpckg(ml, 
                maxbound=riv_rec0.shape[0],
                stress_period_data = riv_spd,
                boundnames = True,
                save_flows=True,
                print_flows=False,
                print_input=False,
                filename = ml.name + '.riv',
                pname = 'riv')

# --- ghb package 
ghb = ml.get_package('ghb_0')
ghb_rec0 = ghb.stress_period_data.get_data()[0]

# fluctuations to reference level

ref_ades_level 
ghb_dh = ades_df.F - ades_df.loc[start_date,'F'] # start_date at low flow, 2023-10-15 (see report)
ghb_dh = ghb_dh.reindex(sim_dates) # re-index with simulation dates
ghb_dh = ghb_dh.interpolate() # gap filling by linear interpolation over simulation period

# gen stress period data 
ghb_spd = {}
for i in range(nper):
    ghb_rec = ghb_rec0.copy()
    ghb_rec['bhead'] = ghb_rec0['bhead']+ghb_dh[i]
    ghb_spd[i]=ghb_rec


ml.ghb.stress_period_data.set_data(ghb_spd)

# --- oc package 
oc = ml.get_package('oc')

spd = { k:[] for k in range(0,nper)}
#for k in range(0,nper,10): spd[k]=[('HEAD','LAST'), ('BUDGET','LAST')] 

oc.saverecord = spd 

# --- IMS package
# loosen convergence criteria to speed up simulation
sim.ims.outer_dvclose = 1e-3
sim.ims.inner_dvclose = 1e-3
sim.ims.print_option='summary'

# ---- ic package 

# load steady state solution
ss_head = flopy.utils.HeadFile(os.path.join('ml','sescousse.hds'))
ss_hdata = ss_head.get_alldata()[0]

# replace initial condition with pre-computed steady state
ml.ic.strt.set_data(ss_hdata)  



# -- write simulation
sim.set_sim_path(transient_dir)
sim.write_simulation()

# -- run simulation 
success, buff = sim.run_simulation(report=True)


# --------------------------------------
# ---- plot input data     ----------
# --------------------------------------

fig,axs = plt.subplots(3,1,sharex=True, figsize=(dbcol_width,dbcol_width)) 

ax0, ax1, ax2 = axs

ax0.bar(p.index, p.values, color='darkblue', label='P')
ax0.bar(et0.index, et0.values, color='tan', label='PET0',alpha=0.8)
ax0.set_ylabel('mm/d')
ax0.legend(loc='upper left')

for pid in ['FS1','FS2','FS3','FS4']:
    obs_levels[pid].plot(ax=ax1,label=pid)

ax1.set_ylabel('m NGF')
ax1.legend(loc='upper right')


for pid in ['PS1','PS2','PS3','FS4']:
    obs_levels[pid].plot(ax=ax2,label=pid)

ax2.set_ylabel('m NGF')
ax2.legend(loc='upper right')


ax2.legend(loc='lower right')

for ax in axs:
    ax.axvline(start_date,ls='--',color='grey')
    ax.axvline(end_date,ls='--',color='grey')

fig.align_ylabels()
fig.tight_layout()

ax.set_xlim(start_date,end_date)

fig.savefig(os.path.join('fig','obs_vars.pdf'),dpi=300)

# --------------------------------------
# ---- plot levels along cross-section      ----------
# --------------------------------------

fig,axs = plt.subplots(3,1,sharex=True, figsize=(10,6)) 

ax0, ax1, ax2, = axs

# weather 
ax0.bar(p.index,p.values, color='darkblue', label='P')
ax0.bar(et0.index,et0.values, color='tan', label='ET0',alpha=0.8)
ax0.set_ylabel('mm/d')
ax0.legend(loc='upper left')

# 1st x-section
pids = ['FS1','FS2','PS1']
linestyles = ['--','--','-']
cols = ['blue','orange','green']
for pid,ls in zip(pids,linestyles):
    l = obs_levels[pid].plot(ax=ax1,ls=ls,label=pid)

ax1.set_ylabel('m NGF')
ax1.legend(loc='upper right')

# 2nd x-section
pids = ['FS2','FS3','PS3']
linestyles = ['--','--','-']
cols = ['blue','orange','green']

for pid,ls in zip(pids,linestyles):
    obs_levels[pid].plot(ax=ax2, ls=ls, label=pid)

ax2.set_ylabel('m NGF')
ax2.legend(loc='upper right')


for ax in axs:
    ax.axvline(start_date,ls='--',color='grey')
    ax.axvline(end_date,ls='--',color='grey')

ax.set_xlim(start_date,end_date)

fig.align_ylabels()
fig.tight_layout()

fig.savefig(os.path.join('fig','piezo_fs.pdf'),dpi=300)

# -----------------------------------------------------------
# ---- investigate regional hydraulic gradient     ----------
# -----------------------------------------------------------


bheads_df = pd.merge( obs_levels[['PS1','PS2','PS3']] ,  ades_df['F'], 
                     left_index=True,right_index=True, 
                     how='left')

fig,axs = plt.subplots(3,1,sharex=True, figsize=(dbcol_width,dbcol_width))
ax0, ax1, ax2 = axs

bheads_df[['PS1','PS2','PS3','F']].plot(ax=ax0)

bheads_df['PS1m'] = bheads_df.PS1 - bheads_df.PS1.mean()
bheads_df['PS2m'] = bheads_df.PS2 - bheads_df.PS2.mean()
bheads_df['PS3m'] = bheads_df.PS3 - bheads_df.PS3.mean()
bheads_df['Fm'] = bheads_df.F - bheads_df.F.mean()

bheads_df[['PS1m','PS2m','PS3m','Fm']].plot(ax=ax1)

dist_PS1_F = 5400 # m 
bheads_df['gradh'] = (bheads_df.PS1 - bheads_df.F) /dist_PS1_F
bheads_df['gradh'].plot(ax=ax2)
ax2.set_ylim(0.0005,0.0007)

fig.align_ylabels()
fig.tight_layout()

fig.savefig(os.path.join('fig','gradh.pdf'),dpi=300)


