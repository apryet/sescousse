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

# load climate forcings 
stjean_file=os.path.join('..','data','stjean_daily.csv')
stjean_df = pd.read_csv(stjean_file, header=[0,1], index_col=0, parse_dates=True)

# ETP from Meteo-France
mf_file=os.path.join('..','data','ETP_daily.csv')
mf_df = pd.read_csv(mf_file,parse_dates=True,header=[0],sep=';')
mf_df.index = pd.to_datetime(mf_df.date,format='%d/%m/%Y')

P = stjean_df.loc[:,('Rain_mm_Tot','sum')]
P.name = 'P'
PET = mf_df.loc[:,'ETP']


# -------------------------------------------------
# model for FS4 
# -------------------------------------------------

# load measured gw levels at FS4
levels_file = os.path.join('..','data','levels_daily.csv')
levels = pd.read_csv(levels_file,header=[0,1],index_col=0,parse_dates=True)
obs_levels = levels.xs('h',axis=1)['FS4'].dropna()


# calibration period
start_date = pd.to_datetime('2023-10-15')
end_date = pd.to_datetime('2024-07-06')

#model setup
ml = ps.Model(obs_levels.loc[start_date:end_date], name="FS4")

sm = ps.RechargeModel(prec=P,evap=PET,
                      rfunc=ps.rfunc.FourParam(), 
                       name="recharge",  
                       recharge=ps.rch.Peterson(),
                       settings=("prec", "evap")
                       ) 

ml.add_stressmodel(sm)

ml.solve(tmin=start_date,tmax=end_date)
fig = ml.plot().get_figure()


fig.savefig(os.path.join('fig','cal_FS4.pdf'),dpi=300)


biais = ml.residuals().mean()
print(biais)

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
start_date = pd.to_datetime('2023-10-15')
end_date = pd.to_datetime('2024-07-06')

#model setup
ml = ps.Model(ades.loc[start_date:end_date], name="ADES")

sm = ps.RechargeModel(prec=P,evap=PET,
                      rfunc=ps.rfunc.Exponential(),
                       name="recharge",  
                       recharge=ps.rch.Linear(),
                       settings=("prec", "evap"),
                       ) 

ml.add_stressmodel(sm)

ml.solve(tmin=start_date,tmax=end_date)
fig = ml.plot().get_figure()


fig.savefig(os.path.join('fig','cal_ADES.pdf'),dpi=300)

biais = ml.residuals().mean()
print(biais)


# -------------------------------------------------
# simulation 
# -------------------------------------------------

mlsim = ps.Model(gw_df.loc[tmin:tmax,'h'], name="GWL")

sm_drias = ps.RechargeModel(prec=drias_df['P'],evap=drias_df['PET'],
                       rfunc=ps.Exponential(), 
                       name="recharge",  
                       recharge=ps.rch.Linear(),
                       settings=("prec", "evap")
                       ) 

mlsim.add_stressmodel(sm_drias)

mlsim.parameters = ml.parameters

sim_gw = mlsim.simulate(tmin='2006-01-01',tmax='2100-12-31')

fig,ax = plt.subplots(1,1,figsize=(6,6))
sim_gw.plot(ax=ax)
gw_df['h'].plot(ls='',marker='.',c='red',ax=ax)

sim_gw.to_csv('gw_RCP85.csv')



