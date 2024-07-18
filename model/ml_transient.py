import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import swb

# --------------------------------------
# ---- settings and data      ----------
# --------------------------------------

# steady model dir
steady_dir = 'ml'

# transient model dir
transient_dir = 'ml_transient'

# load steady state model 
sim = flopy.mf6.MFSimulation.load(sim_ws=steady_dir)
ml = sim.get_model()

# set sim start end dates 
#start_date = pd.to_datetime('2023-07-12')
start_date = pd.to_datetime('2023-10-15')
#end_date = pd.to_datetime('2024-04-12')
end_date = pd.to_datetime('2024-07-06')
sim_dates = pd.date_range(start_date,end_date) 

# load river level recorded at FS4
levels_file = os.path.join('..','data','levels_daily.csv')
levels = pd.read_csv(levels_file,header=[0,1],index_col=0,parse_dates=True)
riv_ref_records = levels.loc[start_date:end_date,('h','FS4')]

# plot to think about the consistency of all that 
df = levels.xs('h',axis=1)
ax = df[['PS1','PS2','PS3','FS4']].plot()
ax.set_xlim(start_date,end_date)
ax.set_ylim(19.5,22.5)
ax.axhline((17.2+21.9)/2,ls='--',color='grey')

# load initial parameter values 
pdata = pd.read_excel(os.path.join('..','data','par.xlsx'),index_col='name')
parvals = pdata['val']

# --------------------------------------
# --- setup recharge model 
# --------------------------------------

# Weather data 
stjean_file=os.path.join('..','data','stjean_daily.csv')
stjean_df = pd.read_csv(stjean_file, header=[0,1], index_col=0, parse_dates=True)

# ETP from Meteo-France
mf_file=os.path.join('..','data','ETP_daily.csv')
mf_df = pd.read_csv(mf_file,parse_dates=True,header=[0],sep=';')
mf_df.index = pd.to_datetime(mf_df.date,format='%d/%m/%Y')

# subset 
P = stjean_df.loc[sim_dates,('Rain_mm_Tot','sum')].values
PET = mf_df.loc[sim_dates,'ETP'].values

# save to single clim file 
clim_df = pd.DataFrame({'P':P,'PET':PET},index=sim_dates)
clim_df.to_csv(os.path.join(transient_dir,'clim.csv'))

swb.run_swb(theta_sat = parvals.loc['tsat'],
            D_max= parvals.loc['dmax'],
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
nper = (end_date - start_date).days
tdis.nper = nper
perlen = 86400 # s
tdis.perioddata = [ [perlen,1,1] for _ in range(nper)]

# --- recharge package
rcha = ml.get_package('rcha_0')
rcha.recharge = {i:rech[i] for i in range(nper)}

# --- evt package (aquifer root water uptake)
evta = flopy.mf6.ModflowGwfevta(
    ml,
    pname="evt",
    save_flows=True,
    surface = 5, # non-limited evt elevation 
    depth = -20, # exctinction elevation 
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
ml.drn.stress_period_data.set_data({0:drn_rec0})

# --- riv package (converted as second drain package drn)
riv = ml.get_package('riv_0')
# original recarray ModflowGwfriv
riv_rec0 = riv.stress_period_data.get_data()[0]

# new ref recarray for ModflowGwfdrn
driv_rec0 = flopy.mf6.ModflowGwfdrn.stress_period_data.empty(ml,maxbound=riv_rec0.shape[0])[0]
driv_rec0['cellid'] = riv_rec0['cellid']
driv_rec0['elev'] = riv_rec0['stage']
driv_rec0['cond'] = parvals.loc['criv']

# reference level
riv_ref_value = 20.6 # min(h_fs4(t))
# fluctuations to reference level
riv_dh = riv_ref_records.values - riv_ref_value

# gen stress period data 
driv_rec = driv_rec0.copy()
driv_spd = {}
for i in range(nper):
    driv_rec['elev'] = driv_rec0['elev']+riv_dh[i]
    driv_spd[i]=driv_rec

driv = flopy.mf6.ModflowGwfdrn(ml, 
                                maxbound=riv_rec0.shape[0],
                                stress_period_data = driv_spd,
                                save_flows=True,
                                print_flows=False,
                                print_input=False,
                                filename = ml.name + '.driv',
                                pname = 'driv')
ml.remove_package('riv')

# --- ghb package 
#ghb = ml.get_package('ghb_0')
#ghb_rec0 = ghb.stress_period_data.get_data()[0]
#ghb_rec0['cond'] = parvals.loc['cghb']
#ml.ghb.stress_period_data.set_data({0:ghb_rec0})

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

# --- save observed levels over simulation period
obs_locs = ['FS1','FS2','FS3','FS4','PS1','PS2','PS3']
obs_levels = levels.xs('h',axis=1)[obs_locs]
obs_levels.index = obs_levels.index.date
obs_levels = obs_levels.loc[start_date:end_date]
obs_levels.to_csv(os.path.join('ml_transient','obs_levels.csv'))

# -- write simulation
sim.set_sim_path(transient_dir)
sim.write_simulation()

# -- run simulation 
success, buff = sim.run_simulation(report=True)
