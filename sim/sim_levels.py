import os
import pandas as pd
import matplotlib.pyplot as plt
import pastas as ps
import numpy as np
from datetime import date, time, datetime

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# -------------------------------------------------
# load data
# -------------------------------------------------

# P from Meteo-France (Mérignac Station)
p_file=os.path.join('..','data','P_daily.csv')
p_df = pd.read_csv(p_file, header=[0],sep=',')
p_df.index = pd.DatetimeIndex(pd.to_datetime(p_df.date,format='%d/%m/%Y'))
P = p_df['P']

# ET0 from Meteo-France (Mérignac Station)
et0_file=os.path.join('..','data','ET0_daily.csv')
et0_df = pd.read_csv(et0_file,header=[0],sep=',')
et0_df.index = pd.DatetimeIndex(pd.to_datetime(et0_df.date,format='%d/%m/%Y'))
ET0 = et0_df['ET0']

# -------------------------------------------------
# model for FS4 
# -------------------------------------------------

# load measured gw levels at FS4
levels_file = os.path.join('..','data','levels_daily.csv')
levels = pd.read_csv(levels_file,header=[0,1],index_col=0,parse_dates=True)
FS4 = levels.xs('h',axis=1)['FS4'].dropna()

# set values to nan when probe is dry
idx = (FS4.index > pd.to_datetime('2024-08-02')) & (FS4.index < pd.to_datetime('2024-09-24'))
FS4.loc[idx] = np.nan

# calibration period
start_date = pd.to_datetime('2023-10-15')
end_date = pd.to_datetime('2024-10-01')

#model setup
ml_riv = ps.Model(FS4.loc[start_date:end_date], name="FS4")

sm_riv = ps.RechargeModel(prec=P,evap=ET0,
                      rfunc=ps.rfunc.FourParam(), 
                       name="recharge",  
                       recharge=pqs.rch.Peterson(),
                       settings=("prec", "evap")
                       ) 

ml_riv.add_stressmodel(sm_riv)

ml_riv.solve(tmin=start_date,tmax=end_date)
fig = ml_riv.plot().get_figure()


fig.savefig(os.path.join('fig','cal_FS4.pdf'),dpi=300)


biais = ml_riv.residuals().mean()
print(biais)

# simulation over  historical period 
sim_FS4 = ml_riv.simulate(tmin='1950-01-01',tmax='2024-10-14')
sim_FS4.name = 'h'
sim_FS4.to_csv('sim_FS4.csv')

# -------------------------------------------------
# model for F ADES 
# -------------------------------------------------

# Piezometric records from ADES
# https://ades.eaufrance.fr/Fiche/PtEau?code=07545X0029/F
ades_file=os.path.join('..','data','BSS001VYWT.csv')
ades_df = pd.read_csv(ades_file,header=[0],sep=',')
ades_df.index = pd.to_datetime(ades_df.date,format='%d/%m/%Y')
ades = ades_df['F'].sort_index()

# calibration period
start_date = pd.to_datetime('2019-01-23')
end_date = pd.to_datetime('2024-10-14')

#model setup
ml_gw = ps.Model(ades.loc[start_date:end_date], name="ADES")

sm = ps.RechargeModel(prec=P,evap=ET0,
                      rfunc=ps.rfunc.Exponential(),
                       name="recharge",  
                       recharge=ps.rch.Linear(),
                       settings=("prec", "evap"),
                       ) 

ml_gw.add_stressmodel(sm)

ml_gw.solve(tmin=start_date,tmax=end_date)
fig = ml_gw.plot().get_figure()

fig.savefig(os.path.join('fig','cal_ADES.pdf'),dpi=300)

biais = ml_gw.residuals().mean()
print(biais)

# simulation over  historical period 
sim_ADES = ml_gw.simulate(tmin='1950-01-01',tmax='2024-10-14')
sim_ADES.name = 'h'
sim_ADES.to_csv('sim_ADES.csv')

# -------------------------------------------------
# simulation over  prospective period 
# -------------------------------------------------

import glob

cm_df = pd.read_excel(os.path.join('..','data','clim_models.xlsx'),index_col=0)
src_dir = os.path.join('..','data','DRIAS')

tmin, tmax = pd.to_datetime('2010-01-01'), pd.to_datetime('2100-01-01')

for cm_id,cm_tag in zip(cm_df.index,cm_df.CM):
    # read processed DRIAS files corresponding to cm_id 
    cm_files = [f for f in glob.iglob(os.path.join(src_dir,f'*{cm_tag}*'))]
    # concatenate historical (<2005) and prospective (>2005) period
    clim = pd.concat([ pd.read_csv(f,sep='\t',parse_dates=True,index_col=0) for f in cm_files])
    clim = clim.sort_index()
    clim.columns = ['P','PET']
    # save clim data
    clim.to_csv(os.path.join('prosp',f'clim_cm{cm_id:02d}.csv'))
    # simulate FS4 with pastas
    sm_riv = ml_riv.stressmodels['recharge']
    sm_riv.prec = ps.timeseries.TimeSeries(df.P)
    sm_riv.evap = ps.timeseries.TimeSeries(df.PET)
    prosp_FS4 = ml_riv.simulate(tmin=tmin,tmax=tmax)
    prosp_FS4.to_csv(os.path.join('prosp',f'sim_FS4_cm{cm_id:02d}.csv'))
    # simulate ADES with pastas 
    sm_gw = ml_gw.stressmodels['recharge']
    sm_gw.prec = ps.timeseries.TimeSeries(df.P)
    sm_gw.evap = ps.timeseries.TimeSeries(df.PET)
    prosp_FS4 = ml_gw.simulate(tmin=tmin,tmax=tmax)
    prosp_FS4.to_csv(os.path.join('prosp',f'sim_ADES_cm{cm_id:02d}.csv'))




