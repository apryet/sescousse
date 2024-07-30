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
gwf.oc.saverecord = [("HEAD", "LAST")]

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
# total drainage 
helpers.write_tot_drn(ws=tpl_dir)

_ = pf.add_observations('sescousse.head.tdiff.csv', 
                            insfile="heads.tdiff.csv.ins", 
                            index_cols="time", 
                            prefix="hdsfluct") 

# add cumulative flows from swb 
_ = pf.add_observations('swb_cum.csv', 
                            insfile="swb_cum.csv.ins", 
                            index_cols="time", 
                            prefix='fcst') 

# add total drainage flow  
_ = pf.add_observations('tot_drn.csv', 
                            insfile="tot_drn.csv.ins", 
                            index_cols="time", 
                            prefix='fcst') 

# add model command
if len(pf.mod_sys_cmds)==0:
    pf.mod_sys_cmds.append('mf6')

# add pre and post- proc functions to forward run 
pf.add_py_function('swb.py', 'run_swb(par_file="par.dat")', is_pre_cmd=True)
pf.add_py_function('helpers.py', 'set_mf_par_vals()', is_pre_cmd=True)
pf.add_py_function('helpers.py', 'process_secondary_obs(ws=".")', is_pre_cmd=False)
pf.add_py_function('helpers.py', 'write_tot_drn()', is_pre_cmd=False)

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

# adjust parameter groups, transformation, bounds and values
par = pst.parameter_data
par.loc[par.index,'pargp'] = pdata.loc[par.index,'pargp']
par.loc[par.index,'partrans'] = pdata.loc[par.index,'partrans']
par.loc[par.index,'parlbnd'] = pdata.loc[par.index,'priorlbnd']
par.loc[par.index,'parubnd'] = pdata.loc[par.index,'priorubnd']
par.loc[par.index,'parval1'] = pdata.loc[par.index,'val']

# update parameter groups
pst.rectify_pgroups()

# generate and save prior parameter cov matrix
pcov = pyemu.Cov.from_parameter_data(pst)
pcov.to_ascii(os.path.join(tpl_dir,'prior_par.cov'))

# generate parmeter realizations
prior_pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=pcov, num_reals=40)
prior_pe.enforce()
prior_pe.to_binary(os.path.join(tpl_dir,'prior_pe.jcb'))

# loosen parameter bounds
par.loc[par.index,'parlbnd'] = pdata.loc[par.index,'parlbnd']
par.loc[par.index,'parubnd'] = pdata.loc[par.index,'parubnd']

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

# fix issue with fs3
idx = obs.index.str.contains('fs3') & (obs.date < pd.to_datetime('2023-12-01').date())
obs.loc[idx,'weight']=0

# less weight to drain obs, anyway
idx = obs.index.str.contains('fs')
obs.loc[idx,'weight'] = 0 # obs.loc[idx,'weight'].div(10)

# 0-weight for forecasts
idx = obs.index.str.contains('fcst')
obs.loc[idx,'weight'] = 0

# fix issue with ps3 (keep only fluctuations)
#idx = obs.index.str.contains('hds_otype:lst_usecol:ps3')
#obs.loc[idx,'weight']=0

# --- further PEST settings 

pst.pestpp_options['uncertainty']='False'

# regularization settings (for regularization mode)
#pst.reg_data.phimlim = pst.nnz_obs
pst.reg_data.phimlim = 3e4 
pst.reg_data.phimaccept = pst.reg_data.phimlim*1.1
pst.reg_data.fracphim = 0.1
pst.reg_data.wfmin = 1.0e-10
pst.reg_data.wfinit = 1e-3
pst.reg_data.wfac = 1.5
pst.reg_data.wtol = 1.0e-2

# SVD parameters 
pst.svd_data.maxsing = pst.npar_adj
pst.svd_data.eigthresh = 1e-6

# adjust derinc for jacobian computation
pst.parameter_groups.loc[:,'derinc']=0.10
pst.parameter_groups.loc[:,'forcen']='always_3'

# implement prefered value regularization
pyemu.helpers.zero_order_tikhonov(pst)

# set forecasts 
forecasts = obs.index[obs.index.str.contains('fcst')].values
pst.pestpp_options['forecasts'] = ','.join(forecasts)

# activate the regularized-GLM solution
pst.pestpp_options['glm_normal_form'] = 'prior'

# define prior parmeter covariance matrix
pst.pestpp_options['parcov'] = 'prior_par.cov'

# set mode to estimation for PESTPP regularized-GLM
pst.control_data.pestmode= 'estimation'

# IES settings 
pst.pestpp_options['ies_num_reals'] = prior_pe.shape[0]
pst.pestpp_options['ies_parameter_ensemble'] = 'prior_pe.jcb'


# set noptmax
pst.control_data.noptmax=0

# write pst 
pst.write(os.path.join(tpl_dir,'cal.pst'))


