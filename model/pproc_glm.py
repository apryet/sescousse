import os, sys, glob
import shutil
import pandas as pd
import numpy as np
import pyemu
import matplotlib.pyplot as plt
import matplotlib
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

# ---------------------------------------------
# load pst
# ---------------------------------------------
phi = pd.read_csv(os.path.join('cal.iobj'),index_col=0)

fig,ax = plt.subplots(1,1,figsize=(sgcol_width,sgcol_width))
phi['measurement_phi'].plot(legend=False,lw=2,color='k', ax=ax,label='Measurement $\Phi$')
ax.set_xlabel('Itérations GLM')
ax.set_ylabel(r'Measurement $\Phi$')
lhm,llm = ax.get_legend_handles_labels(
twax = ax.twinx()
phi['regularization_phi'].plot(legend=False,lw=2,ls='--',color='grey', ax=twax,label='Regul. $\Phi$')
twax.set_ylabel(r'Regularization $\Phi$')
lhr,llr = twax.get_legend_handles_labels()

ax.legend(lhm+lhr,llm+llr,loc='upper right')

fig.tight_layout()

fig.savefig(os.path.join('fig','phi_evol.pdf'),dpi=300)

# ---------------------------------------------
# parameters
# ---------------------------------------------

parnmes=['k','sy','cdrn','criv', 'cghbe', 'cghbw', 'dmax', 'tsat','kc']
#parnmes = pst.parameter_data.index

parsum  = pd.read_csv('cal.par.usum.csv',index_col=0)

#  --- plot sim-obs scatter plots 
fig,axs = plt.subplots(3,4,figsize=(mdcol_width,mdcol_width))

for ax,parnme in zip(axs.ravel(),parnmes) :
    pr_val = np.log10(pst.parameter_data.loc[parnme,'parval1'])
    pr_le = parsum.loc[parnme,'prior_mean'] - parsum.loc[parnme,'prior_lower_bound']
    pr_ue = parsum.loc[parnme,'prior_upper_bound']- parsum.loc[parnme,'prior_mean']
    pt_val = parsum.loc[parnme,'post_mean']
    pt_le = parsum.loc[parnme,'post_mean'] - parsum.loc[parnme,'post_lower_bound']
    pt_ue = parsum.loc[parnme,'post_upper_bound']- parsum.loc[parnme,'post_mean']
    ax.errorbar([0,1], [pr_val,pt_val], 
                yerr=[[pr_le,pt_le],[pr_ue,pt_ue]],
                markersize=6,capsize=4,
                elinewidth = 2,
                color='r',ecolor='k',fmt='.',label='prior')
    ax.set_xticks([0,1],['pr','pt'])
    ax.set_xlim(-1,2)
    ax.set_title(parnme)

fig.tight_layout()
fig.savefig(os.path.join('fig','par_usum.pdf'),dpi=300)

# ---------------------------------------------
# predictions
# ---------------------------------------------
predsum = pd.read_csv('cal.pred.usum.csv')
predsum.index= ['RU','R','DRN','T']
predsum['pr_le'] = predsum['prior_mean'] - predsum['prior_lower_bound']
predsum['pr_ue'] = predsum['prior_upper_bound'] - predsum['prior_mean']
predsum['pt_le'] = predsum['post_mean'] - predsum['post_lower_bound']
predsum['pt_ue'] = predsum['post_upper_bound'] - predsum['post_mean']

fig,ax = plt.subplots(1,1,figsize=(sgcol_width,sgcol_width))
predsum['prior_mean'].plot(ax=ax,kind='bar',yerr=[predsum['pt_le'],predsum['pt_ue']],capsize=4)
ax.set_ylabel('Hauteur d\'eau cumulée [mm]')
fig.tight_layout()
fig.savefig(os.path.join('fig','pred_usum.pdf'),dpi=300)

