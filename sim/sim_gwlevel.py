import os
import pandas as pd
import matplotlib.pyplot as plt
import pastas as ps
import numpy as np
from datetime import date, time, datetime


# -------------------------------------------------
# load data
# -------------------------------------------------

# load measured weather at Mérignac
metfile = os.path.join('data','merignac4.xlsx')
met_df = pd.read_excel(metfile)
met_df.set_index(pd.to_datetime(met_df.DATE,format='%Y%m%d'),inplace=True)
met_df['P']=met_df['RR']
met_df['PET']=met_df['ETPMON']
# load measured gw levels at Hourtin
gwlevel_file = os.path.join('data','BSS001VYWT.csv')
gw_df = pd.read_csv(gwlevel_file)
gw_df.set_index(pd.to_datetime(gw_df.date,format='%d/%m/%Y'),inplace=True)

# load DRIAS scenario

# reference (1951-2005)

ref_file = os.path.join('data','P07327_REF.txt')
ref_df = pd.read_csv(ref_file,header=None,skiprows=96,
                      index_col=0,usecols=[0,3,4,10],
                      parse_dates=True)

# RCP8.5  (2006-2100)
rcp85_file = os.path.join('data','P07327_RCP85.txt')
rcp85_df = pd.read_csv(rcp85_file,header=None,skiprows=96,
                      index_col=0,usecols=[0,3,4,10],
                      parse_dates=True)

drias_df = pd.concat((ref_df,rcp85_df),axis=0)
drias_df.index.name = 'date'
drias_df.columns = ['T','P','PET']
drias_df['PET'] = drias_df.PET*86400

drias_df.to_excel('comp_drias.xlsx')

met_df[['P','PET']].plot()
drias_df.plot()

# -------------------------------------------------
# ref drias vs met analysis 
# -------------------------------------------------

dmin = pd.to_datetime('1986-01-01')
dmax = pd.to_datetime('2005-12-31')
dates = pd.date_range(dmin,dmax)
drias_ss = drias_df.loc[dates][['P','PET']]
met_ss = met_df.loc[dates][['P','PET']]

fig, axs = plt.subplots(1,2,figsize=(6,4))
axs[0].plot(drias_ss.P,met_ss.P,'+',c='blue')
axs[0].set_aspect('equal')
axs[1].plot(drias_ss.PET,met_ss.PET,'+',c='orange')
axs[1].set_aspect('equal')


# -------------------------------------------------
# model 
# -------------------------------------------------

# calibration period
tmin='2013-10-01'
tmax='2022-03-01'

#model setup
ml = ps.Model(gw_df.loc[tmin:tmax,'h'], name="GWL")

sm = ps.RechargeModel(prec=met_df['P'],evap=met_df['PET'],
                       rfunc=ps.Exponential(), 
                       name="recharge",  
                       recharge=ps.rch.Linear(),
                       settings=("prec", "evap")
                       ) 

ml.add_stressmodel(sm)

ml.solve(tmin=tmin,tmax=tmax)

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


# -------------------------------------------------
# débit matouse 
# -------------------------------------------------

# stream discharge Matouse
qmat_file = os.path.join('data','debit_matouse.csv')
qmat_df = pd.read_csv(qmat_file)
qmat_df.set_index(pd.to_datetime(qmat_df.date,format='%d/%m/%Y'),inplace=True)


# calibration period
tmin='2000-01-01'
tmax='2020-12-31'


ps.Gamma()

Linear, FlexModel and Berendrecht
#model setup
qml = ps.Model(qmat_df.loc[tmin:tmax,'Q'], name="Q")

qsm = ps.RechargeModel(prec=met_df['P'],evap=met_df['PET'],
                       rfunc=ps.Gamma(), 
                       name="recharge",  
                       recharge=ps.rch.Berendrecht(),
                       settings=("prec", "evap")
                       ) 

qml.add_stressmodel(qsm)

qml.solve(tmin=tmin,tmax=tmax)
qml.plot()

