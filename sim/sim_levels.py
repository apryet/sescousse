import os
import pandas as pd
import pastas as ps
import numpy as np
from datetime import date, time, datetime
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob


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
                       recharge=ps.rch.Peterson(),
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

rcp_start= pd.to_datetime('2005-08-01') # rcp start 
cm_df = pd.read_excel(os.path.join('..','data','clim_models.xlsx'),index_col=0)
src_dir = os.path.join('..','data','DRIAS')

tmin, tmax = pd.to_datetime('1980-01-01'), pd.to_datetime('2099-12-31')

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
    sm_riv.prec = ps.timeseries.TimeSeries(clim.P)
    sm_riv.evap = ps.timeseries.TimeSeries(clim.PET)
    prosp_FS4 = ml_riv.simulate(tmin=tmin,tmax=tmax,warmup=5*365)
    prosp_FS4.to_csv(os.path.join('prosp',f'sim_FS4_cm{cm_id:02d}.csv'))
    # simulate ADES with pastas 
    sm_gw = ml_gw.stressmodels['recharge']
    sm_gw.prec = ps.timeseries.TimeSeries(clim.P)
    sm_gw.evap = ps.timeseries.TimeSeries(clim.PET)
    prosp_FS4 = ml_gw.simulate(tmin=tmin,tmax=tmax,warmup=5*365)
    prosp_FS4.to_csv(os.path.join('prosp',f'sim_ADES_cm{cm_id:02d}.csv'))


# -------------------------------------------------
# -----            plots        -------------------
# -------------------------------------------------

prosp_dic = {}

for cm_id in cm_df.index:
    clim_file = os.path.join('prosp',f'clim_cm{cm_id:02d}.csv')
    clim = pd.read_csv(clim_file,header=0,index_col=0,parse_dates=True)
    ades_file = os.path.join('prosp',f'sim_ADES_cm{cm_id:02d}.csv')
    ades = pd.read_csv(ades_file,header=0,index_col=0,parse_dates=True)
    ades.columns=['ADES']
    fs4_file = os.path.join('prosp',f'sim_FS4_cm{cm_id:02d}.csv')
    fs4 = pd.read_csv(fs4_file,header=0,index_col=0,parse_dates=True)
    fs4.columns=['FS4'] 
    prosp_dic[cm_id] = pd.concat([clim,ades,fs4],axis=1)[tmin:tmax]

prosp = pd.concat(prosp_dic,names=['cm'],axis=1)
agg_dic = {'P':'sum','PET':'sum','ADES':'mean','FS4':'mean'}
mix_agg_dic = {(l,c):f for l in prosp.columns.levels[0] for c,f in agg_dic.items()}
prospy = prosp.groupby(pd.Grouper(freq='Y')).agg(mix_agg_dic)


fig,axs=plt.subplots(4,1,sharex=True,figsize=(10,7))
# ---- total precip
tsy=prospy.xs('P',1,1)
tsy.plot(style='+',ax=axs[0],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
tsy.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[0],color='tomato',ls='-', lw=2,legend=False)
axs[0].grid(which='both')
axs[0].set_xticklabels([])
axs[0].set_ylabel('P [mm/y]')

# ---- potential evapotranspiration
tsy=prospy.xs('PET',1,1)
tsy.plot(style='+',ax=axs[1],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
tsy.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[1],color='tomato',ls='-', lw=2,legend=False)
axs[1].grid(which='both')
axs[1].set_xticklabels([])
axs[1].set_ylabel('PET [mm/y]')

# ---- ADES Groundwater level 
tsy=prospy.xs('ADES',1,1)
tsy.plot(style='+',ax=axs[2],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
tsy.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[2],color='tomato',ls='-', lw=2,legend=False)
axs[2].grid(which='both')
axs[2].set_xticklabels([])
axs[2].set_ylabel('ADES [m NGF]')

# ---- FS4 river level 
tsy=prospy.xs('FS4',1,1)
tsy.plot(style='+',ax=axs[3],color='black',lw=0.5,ms=3, alpha=0.5,legend=False)
tsy.mean(axis=1).rolling(window=10,center=True).mean().plot(ax=axs[3],color='tomato',ls='-', lw=2,legend=False)
axs[3].grid(which='both')
axs[3].set_ylabel('FS4 [m NGF]')

# --- decoration
for ax in axs:
    ax.axvline(pd.to_datetime(rcp_start),alpha=0.5,lw=1.5,ls=':',color='black')

lls =  [ Line2D([0], [0], label=f'Simulated Annual Values', marker='+', linestyle='', color='black',alpha=0.5)]
lls += [Line2D([0], [0], label='Multi-model 10-year moving average',linestyle='-', color='tomato')]
fig.legend(handles=lls,loc='upper center',ncols=6,facecolor='white', framealpha=1)
fig.tight_layout()
fig.savefig(os.path.join('fig','long_term_records.pdf'),dpi=300)


