import os, sys, glob
import shutil
import pandas as pd
import numpy as np
import pyemu
import matplotlib.pyplot as plt
import matplotlib
import flopy 
from matplotlib.lines import Line2D

# plot settings
plt.rc('font', family='serif', size=9)
sgcol_width = 9/2.54
mdcol_width = 14/2.54
dbcol_width = 19/2.54

# ---------------------------------------------
# load pst
# ---------------------------------------------

pstfile = 'cal.pst'
case = pstfile.split('.')[0]
pst = pyemu.Pst(pstfile)
forecasts = pst.pestpp_options['forecasts'].split(',')

# iteration id for posterior distribution 
pt_id = pst.control_data.noptmax
pt_id = 3

# load prior and last iteration observation ensembles  
pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.0.obs.csv"))
pr_oe._df.index=pr_oe.index.astype(str) #trick since index mixed types in oe.index (should find a fix)

# posterior (last iteration) ensemble of simulated values 
pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.{pt_id}.obs.csv"))
pt_oe._df.index=pt_oe.index.astype(str) #trick since index mixed types in oe.index (should find a fix)

# observation with noise 
obswns = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(f"{case}.obs+noise.csv"))

# ---------------------------------------------
# load model
# ---------------------------------------------

sim_dir = '.'
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir)
start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())

# ---------------------------------------------
# load prior and last iteration parameter ensembles  
# ---------------------------------------------

# load parameter ensembles 
pe_dic = { it:pyemu.ParameterEnsemble.from_csv(
    pst=pst,filename=os.path.join(f"{case}.{it}.par.csv")
    ) for it in range(pt_id+1)}


# load prior parameter ensemble 
pr_pe = pe_dic[0]

# load posterior parameter ensemble 
pt_pe = pe_dic[pt_id]

# ---------------------------------------------
# evolution of phi 
# ---------------------------------------------
fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10,3.5))
# left
ax = axes[0]
phi = pd.read_csv(os.path.join(f"{case}.phi.actual.csv"),index_col=0)
phi.index = phi.total_runs
phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.5,color='k', ax=ax)
ax.axhline(np.log10(32787.8),ls='--',color='darkgreen',label='GLM')
ax.set_title(r'Actual ')
ax.set_ylabel(r'log ')
# right
ax = axes[-1]
phi = pd.read_csv(os.path.join(f"{case}.phi.meas.csv"),index_col=0)
phi.index = phi.total_runs
phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.2,color='r', ax=ax)
ax.set_title(r'Measured+Noise ')
fig.tight_layout()
fig.savefig(os.path.join('fig','phi_evol.pdf'),dpi=300)

# ---------------------------------------------
# parameters pt - pr distributions 
# ---------------------------------------------

# get parameter df and append phi value value 
pr_df = pr_pe._df.copy(deep=True)
pt_df = pt_pe._df.copy(deep=True)

# load GLM results, for comparison
parsum  = pd.read_csv(os.path.join('..','master_glm','cal.par.usum.csv'),index_col=0)

# plot hydraulic properties (per layer)
parnmes=['k','sy','cdrn','criv', 'cghbe', 'cghbw', 'dmax', 'tsat','kc']
parnmes = pst.parameter_data.index

fig, axs = plt.subplots(3, 3, figsize=(dbcol_width,dbcol_width))

for p,ax in zip(parnmes,axs[:len(parnmes)].ravel()):
    vmin,vmax = np.round(np.log10([pr_df[p].min(),pr_df[p].max()]))
    logbins = np.linspace(vmin-1,vmax+1,10)
    pr_df[p].apply(np.log10).hist(ax=ax,fc='0.5',ec="none",alpha=0.5,
                                         density=False,
                                         bins=logbins,label='prior')
    pt_df[p].apply(np.log10).hist(ax=ax,fc='b',ec="none",alpha=0.5,
                                         density=False,
                                         bins=logbins,label='posterior')
    ax.axvline(parsum.loc[p,'post_mean'],c='red',ls='--')
    ax.set_title(p)

fig.tight_layout()

fig.savefig(os.path.join('fig',f'par_pr_pt_hist.pdf'),dpi=300)

# ---------------------------------------------
# fit to observations 
# ---------------------------------------------

# get ensemble of time series for given observation group 
def get_og_ts(oe,onames,odates, trans):
    ts = oe._df.loc[:,onames].T.apply(trans)
    ts.index = odates 
    return(ts)

def plot_tseries_ensembles(pr_oe, pt_oe, obswns, ognmes, ogdates, trans=None, ylabel='',legend=True ):
    # get the observation data from the control file and select 
    obs = pst.observation_data.copy()
    fig,axes = plt.subplots(len(ognmes),1,sharex=True,figsize=(dbcol_width,0.8*dbcol_width))
    if trans==None:
        trans=[lambda x : x]*len(ognmes)
    if not type(trans)==list:
        trans = [trans]*len(ognmes)
    # for each observation group (i.e. timeseries)
    for ax,og,t in zip(axes,ognmes,trans):
        # get values 
        oobs = obs.loc[obs.obgnme==og.lower(),:].copy()
        onames = oobs.obsnme.values
        odates = ogdates[og]
        # plot prior
        if pr_oe is not None :
            ts = get_og_ts(pr_oe,onames,odates, trans=t)
            ts.plot(ax=ax,color='grey',lw=0.5,alpha=0.40,legend=False)
        # plot posterior
        if pt_oe is not None :
            ts = get_og_ts(pt_oe,onames,odates,trans=t)
            ts.plot(ax=ax,color='red',lw=0.5,alpha=0.40,legend=False)
            ts['base'].plot(ax=ax,color='green',alpha=1,lw=1,legend=False)
        # plot measured+noise 
        if obswns is not None :
            ts = get_og_ts(obswns,onames,odates,trans=t)
            ts.plot(ax=ax,color='blue',lw=0.5,alpha=0.40,legend=False)
        # plot obs
        ax.plot(odates, oobs.obsval.apply(t).values,'+',color="black",ms=3,alpha=0.5,lw=1)
        ax.set_title(og,loc="left")
        ax.set_ylabel(ylabel)
        lpr = Line2D([0], [0], label='Sim. prior', color='grey')
        lpt = Line2D([0], [0], label='Sim. posterior', color='red')
        lbase = Line2D([0], [0], label='Sim. base', color='green')
        lobs = Line2D([0], [0],marker='+',ls='None', label='Observed', color='black')
        lobsn = Line2D([0], [0], label='Obs.+noise', color='blue')
        if legend:
            ax.legend(handles=[lpr,lpt,lbase,lobs,lobsn],loc='upper left',ncols=5)
            plot_legend=False
    fig.tight_layout()
    return fig

# 

obgnmes = ['oname:hds_otype:lst_usecol:ps1','oname:hds_otype:lst_usecol:ps2','oname:hds_otype:lst_usecol:ps3']

# get date sequence per observation groups 
ogdates = {}
for obgnme in obgnmes:
    onmes = pst.observation_data.index[pst.observation_data.obgnme == obgnme]
    otimes = [float(obnme.split('time:')[1]) for obnme in onmes]
    dates = pd.to_datetime( start_date + pd.to_timedelta(otimes,'s'))
    ogdates[obgnme]=dates


# observation wells
fig = plot_tseries_ensembles(pr_oe, pt_oe, obswns , obgnmes, ogdates,trans=lambda x : x, ylabel='Groundwater level [m NGF]')

# reset axis titles 
axs = fig.get_axes()
for ax,obgnme in zip(axs,obgnmes) :
    ax.set_title(obgnme[-3:].upper(),loc='left')


fig.savefig(os.path.join('fig','pr_pt_qsimobs.png'),dpi=300)



# ---------------------------------------------
# forecasts pt - pr distributions 
# ---------------------------------------------


pr_oe_df = pr_oe._df.copy(deep=True)
pt_oe_df = pt_oe._df.copy(deep=True)

pr_df = pr_oe_df.loc[:,forecasts]
pt_df = pt_oe_df.loc[:,forecasts]

axs = pr_df.hist(fc="0.5",ec="none",alpha=0.5,density=False,label='prior')
axs = pt_df.hist(ax=axs.ravel()[:len(parnmes)],fc="b",ec="none",alpha=0.5,density=False,label='posterior')


# reset axis titles 
for ax in axs :
    nme = ax.get_title()
    nme = nme.split('usecol:')[-1].split('_time')[0].upper()
    ax.set_title(nme)

fig = axs.ravel()[0].get_figure()

fig.tight_layout()
fig.savefig(os.path.join('fig',f'forecasts_pr_pt_hist.pdf'),dpi=300)

