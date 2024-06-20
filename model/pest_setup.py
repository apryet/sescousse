import os, shutil
import numpy as np
import pandas as pd
import flopy
import pyemu
import helpers


# -----------------------------------------------------------
# admin 
# -----------------------------------------------------------

# folder containing original model files
org_dir = os.path.join('ml_transient')

# a dir to hold a copy of the org model files
tmp_dir = os.path.join('ml_tmp')

if os.path.exists(tmp_dir):
    shutil.rmtree(tmp_dir)

shutil.copytree(org_dir,tmp_dir)

# load simulation
sim = flopy.mf6.MFSimulation.load(sim_ws=tmp_dir)

# load flow model
gwf = sim.get_model()

# only save last head heads 
gwf.oc.saverecord = [] # [("HEAD", "LAST")]

# write oc package (!)
gwf.oc.write()

# run model
#success, buff = sim.run_simulation(report=True)

# clear heavy output files 
rm_files = [f'{gwf.name}.hds',f'{gwf.name}.list']
for f in rm_files:
    fpath = os.path.join(tmp_dir,f)
    if os.path.exists(fpath):
        os.remove(fpath)

# spatial reference 
sr = pyemu.helpers.SpatialReference.from_namfile(
        os.path.join(tmp_dir, f'{gwf.name}.nam'),
        delr=gwf.dis.delr.array, delc=gwf.dis.delc.array)

# model template dir for pest 
tpl_dir = os.path.join('ml_tpl')
start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())

# instantiate PstFrom
pf = pyemu.utils.PstFrom(original_d=tmp_dir, 
                            new_d=tpl_dir, 
                            remove_existing=True, 
                            longnames=True, 
                            spatial_reference=sr, 
                            zero_based=False, 
                            start_datetime=start_date, 
                            echo=False)
     
# -- add observations 

# absolute head values 
hsim_file = 'sescousse.head.csv'
_ = pf.add_observations(hsim_file, 
                            insfile=f'{hsim_file}.ins', 
                            index_cols='time', 
                            prefix='hds')
# head fluctuations to mean
helpers.process_secondary_obs(ws=tpl_dir)

_ = pf.add_observations('sescousse.head.tdiff.csv', 
                            insfile="heads.tdiff.csv.ins", 
                            index_cols="time", 
                            prefix="hdsfluct") 

# add model command
if len(pf.mod_sys_cmds)==0:
    pf.mod_sys_cmds.append('mf6')

# add pre and post- proc functions to forward run 
pf.add_py_function('swb.py', 'run_swb(par_file="par.dat")', is_pre_cmd=True)
pf.add_py_function('helpers.py', 'set_mf_par_vals()', is_pre_cmd=True)
pf.add_py_function('helpers.py', 'process_secondary_obs(ws=".")', is_pre_cmd=False)

# -- generate pst 
pst = pf.build_pst()

# --- process parameters

# load parameter settings from excel 
pdata = pd.read_excel(os.path.join('..','data','par.xlsx'),index_col='name')

par_file = os.path.join(tpl_dir,'par.dat')
tpl_file = par_file + '.tpl'

#  write template file 
with open(tpl_file, "w") as f:
    parnames = pdata.index.values
    f.write("ptf ~\n")
    [f.write("{0:<12}~{0:^12}~\n".format(cname)) for cname in parnames]


# add parameters 
par = pst.add_parameters(tpl_file,pst_path='.')

# adjust parameter bounds and values
par = pst.parameter_data
par.loc[par.index,'parlbnd'] = pdata.loc[par.index,'parlbnd']
par.loc[par.index,'parubnd'] = pdata.loc[par.index,'parubnd']
par.loc[par.index,'parval1'] = pdata.loc[par.index,'val']

# --- process observed values and weights 

# load sim values
hsim = pd.read_csv(os.path.join(tpl_dir,'sescousse.head.csv'),index_col='time')
locs = hsim.columns
hsim['date'] = (pd.to_datetime(start_date)+pd.to_timedelta(hsim.index.values.astype(float),'s')).date

# load observed values in a single, multiindexed series 
levels_file = os.path.join('..','data','levels_daily.csv')
levels = pd.read_csv(levels_file, header=[0,1], index_col=0, parse_dates=True)
hobs = levels.xs('h',axis=1)[locs]
hobs.index = hobs.index.date # convert datetime to date 
hobs = hobs.stack() # convert to multindexed series of heads
hobs.index.names = ['date','loc']
hobs = pd.concat([hobs],keys=['hds'],names=['oname'])

# compute observed head fluctuations and convert to a single, multiindexed series
hobsfluct = levels.xs('h',axis=1)[locs] - levels.xs('h',axis=1)[locs].mean()
hobsfluct.index = hobsfluct.index.date # convert datetime to date 
hobsfluct = hobsfluct.stack() # convert to multindexed series of heads
hobsfluct.index.names = ['date','loc']
hobsfluct = pd.concat([hobsfluct],keys=['hdsfluct'],names=['oname'])

# expected sim obs misfit 
sigma = 0.10 # m

# add date columns to pst.observation_data 
obs = pst.observation_data
obs['date'] = (pd.to_datetime(start_date)+pd.to_timedelta(obs.time.values.astype(float),'s')).date
# temporarily reset obs index 
obs.index = pd.MultiIndex.from_frame(pd.DataFrame({'oname':obs.oname,'date':obs.date, 'loc':obs.usecol.str.upper()}))
# default weights to 0 to make sure unavailable obs are discarded in PHI
obs.weight = 0
# set head obs 
idx = obs.index.intersection(hobs.index) 
obs.loc[idx,'obsval'] = hobs.loc[idx]
obs.loc[idx,'weight'] = 1/sigma
# set head fluct obs
idx = obs.index.intersection(hobsfluct.index) 
obs.loc[idx,'obsval'] = hobsfluct.loc[idx]
obs.loc[idx,'weight'] = 1/sigma

# reset index 
obs.index= obs.obsnme

# set weight to 0 when drains are dry 
idx = obs.index.str.contains('fs') & (obs.date < pd.to_datetime('2023-11-01').date())
obs.loc[idx,'weight']=0

# --- further PEST settings 

# regularization settings
pst.reg_data.phimlim = pst.nnz_obs
pst.reg_data.phimaccept = pst.reg_data.phimlim*1.1
pst.reg_data.fracphim = 0.1
pst.reg_data.wfmin = 1.0e-10
pst.reg_data.wfinit = 1e-3
pst.reg_data.wfac = 1.5
pst.reg_data.wtol = 1.0e-2

# adjust derinc
pst.parameter_groups.loc[:,'derinc']=0.05

# implement prefered value regularization
pyemu.helpers.zero_order_tikhonov(pst)

# set noptmax
pst.control_data.noptmax=0

# write pst 
pst.write(os.path.join(tpl_dir,'cal.pst'))


