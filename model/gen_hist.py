import sys, os, shutil
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import swb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tpl_dir')
parser.add_argument('--sim_dir')
args = parser.parse_args()

tpl_dir = args.tpl_dir
sim_dir = args.sim_dir

# template dir 
#tpl_dir = 'drn110'

# simulation directory
#sim_dir = 'histo_drn110'

# cp simulation from template 
if not os.path.exists(sim_dir):
    shutil.copytree(tpl_dir,sim_dir)
else: 
    print('***WARNING***: Failed to copy, directory exists')

# set sim start end dates 
start_date = pd.to_datetime('2004-10-01').date()
end_date = pd.to_datetime('2024-09-30').date()
sim_dates = pd.date_range(start_date,end_date).date

# --------------------------------------------------
# load data 
# --------------------------------------------------

# load simulation forcing variables 
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
clim_df.ET0 = clim_df.ET0.interpolate() # gap-filling
clim_df.to_csv(os.path.join(sim_dir,'clim.csv'))

# simulated piezometric level at ADES obs well
ades_file=os.path.join('..','sim','histo','sim_ADES.csv')
ades_df = pd.read_csv(ades_file,parse_dates=True,index_col=0,sep=',')
ades_df.index = pd.DatetimeIndex(pd.to_datetime(ades_df.index,format='%d/%m/%Y')).date
ades = ades_df['h']

# simulated river level at FS4 
fs4_file=os.path.join('..','sim','histo','sim_FS4.csv')
fs4_df = pd.read_csv(fs4_file,parse_dates=True,index_col=0,sep=',')
fs4_df.index = pd.DatetimeIndex(pd.to_datetime(fs4_df.index,format='%d/%m/%Y')).date
fs4 = fs4_df['h']

# --------------------------------------------------
# update stress period data 
# --------------------------------------------------

# load template model 
sim = flopy.mf6.MFSimulation.load(sim_ws=tpl_dir)
ml = sim.get_model()

# --- setup recharge model 

parvals = pd.read_csv(os.path.join(sim_dir,'par.dat'),header=None,delim_whitespace=True,index_col=0)[1]

swb.run_swb(theta_sat = parvals.loc['tsat'],
            D_max= parvals.loc['dmax'],
            Kc = parvals.loc['kc'],
            cwd=sim_dir)

swb_df = pd.read_csv(os.path.join(sim_dir,'swb_vars.csv'),index_col=0,parse_dates=True)
rech = swb_df.loc[:,'R']*0.001/86400 # mm/d to m/s # recharge
evt = swb_df.loc[:,'RU']*0.001/86400 # mm/d to m/s # transpiration (root water uptake)

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
evta = ml.get_package('evt')
evta.rate.set_data({i:evt[i] for i in range(nper)})

# --- drn package (if present)
try : 
    drn = ml.get_package('drn_0')
    drn_spd = drn.stress_period_data.get_data()
    rec0 = drn_spd[0]
    drn.stress_period_data.set_data({i:rec0 for i in range(nper)})
except : # in case the scenario does not present drn package 
    print(f'***INFO***: Could not update stress_period_data for {sim_dir}')

# --- river network 

# river level records at FS4
riv_ref_records = fs4.loc[start_date:end_date]

# load package 
riv = ml.get_package('riv')
riv_rec0 = riv.stress_period_data.get_data(0)

# reference value at FS4 for steady state simulation
# obs_levels.loc[pd.to_datetime('2023-10-15').date(),'FS4']
riv_ref_value = 20.66 # hriv(fs4,ss)

# gen stress period data 
riv_spd = {}
for i in range(nper):
    riv_rec = riv_rec0.copy()
    # h_riv = hriv(ss) + (hriv(fs4,t)-hriv(fs4,ss)) + dh
    riv_rec['stage'] = riv_rec0['stage']+ (riv_ref_records[i]-riv_ref_value)+parvals.loc['dhriv']
    # set conductance 
    riv_rec['cond'] = parvals.loc['criv']
    # adjust river bottom to avoid hriv < rbot
    riv_rec['rbot'] = np.minimum(riv_rec['stage']-0.01,riv_rec0['rbot'])
    riv_spd[i]=riv_rec

riv.stress_period_data.set_data(riv_spd)

# --- ghb package 
ghb = ml.get_package('ghb_0')
ghb_rec0 = ghb.stress_period_data.get_data()[0]

# fluctuations to reference level
ghb_dh = ades_df.h - 16.69  # ades_df.loc[tref,'F']=16.69, with tref
ghb_dh = ghb_dh.loc[start_date:end_date]

# gen stress period data 
ghb_spd = {}
for i in range(nper):
    ghb_rec = ghb_rec0.copy()
    ghb_rec['bhead'] = ghb_rec0['bhead'] + ghb_dh[i]
    ghb_spd[i]=ghb_rec

ml.ghb.stress_period_data.set_data(ghb_spd)

# --- oc package 
oc = ml.get_package('oc')

spd = { k:[] for k in range(0,nper)}
for k in range(0,nper): spd[k]=[('HEAD','LAST')] 

oc.saverecord = spd 

# -- write simulation
sim.set_sim_path(sim_dir)
sim.write_simulation()

# -- run simulation 
success, buff = sim.run_simulation(report=True)


