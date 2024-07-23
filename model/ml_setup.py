import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import flopy
from mfsetup import MF6model

wd = os.getcwd()

# --- load initial parameter values and update parameter values 
pdata = pd.read_excel(os.path.join('..','data','par.xlsx'),index_col='name')
parvals = pdata['val']

# generate model 
m = MF6model.setup_from_yaml('sescousse.yml',verbose=True)
m.write_input()

# write grid to shapefile 
m.modelgrid.write_shapefile('postproc/shps/grid.shp')

# load simulation and run model 
sim = flopy.mf6.MFSimulation.load(wd)
ml = sim.get_model('sescousse')


# set k 
ml.npf.k.set_data(parvals.loc['k'])

# cdrn  
drn = ml.get_package('drn_0')
drn_rec0 = drn.stress_period_data.get_data()[0]
drn_rec0['cond'] = parvals.loc['cdrn']
ml.drn.stress_period_data.set_data({0:drn_rec0})

# set criv
riv = ml.get_package('riv_0')
riv_rec0 = riv.stress_period_data.get_data()[0]
riv_rec = riv_rec0.copy()
riv_rec['cond'] = parvals.loc['criv']

# ghb
print('Setting cghb values...')
ghb_rec0 = ml.ghb.stress_period_data.get_data()[0]
# western bc
idx = ghb_rec0['boundname']=='w'
ghb_rec0['cond'][idx] = parvals.loc['cghbw']
ghb_rec0['bhead'][idx] = parvals.loc['hghbw']
# eastern bc
idx = ghb_rec0['boundname']=='e'
ghb_rec0['cond'][idx] = parvals.loc['cghbe']
ghb_rec0['bhead'][idx] = parvals.loc['hghbe']
# update values 
ml.ghb.stress_period_data.set_data({0:ghb_rec0})

# ims (strengten convergence criteria for steady-state)
sim.ims.outer_dvclose = 1e-5
sim.ims.inner_dvclose = 1e-5

# --- run simulation
sim.write_simulation()
success, buff = sim.run_simulation()

# ---- plot figures 

# make fig dir 
if not os.path.exists('fig'):
    os.mkdir('fig')

xlabel = 'x (estward) m RGF93'
ylabel = 'y (northward) m RGF93'
masked_values = [-9999]

# load head 
head = flopy.utils.HeadFile('sescousse.hds')
hdata = head.get_alldata()[0]

# head cross-section
fig,ax = plt.subplots(1,1,figsize=(4,4))
ax.plot(hdata[0,1:-1, int(ml.dis.ncol.data/2)])
fig.tight_layout()
fig.savefig(os.path.join('fig','xs_head.png'),dpi=300)

# plot head 
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Head")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(hdata)#,vmin=19,vmax=21)
cb = plt.colorbar(quadmesh,label='m NGF', shrink=0.5, ax=ax)
fig.tight_layout()
fig.savefig(os.path.join('fig','heads.png'),dpi=300)

# plot boundary conditions 
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_title("Model Boundary Conditions")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
pmv = flopy.plot.PlotMapView(model=ml, ax=ax)
linecollection = pmv.plot_grid(linewidth=0.1)
quadmesh = pmv.plot_bc('riv')
quadmesh = pmv.plot_bc('ghb')
quadmesh = pmv.plot_bc('drn')
fig.tight_layout()
fig.savefig(os.path.join('fig','bc.png'),dpi=300)

# plot top array
top = ml.dis.top.array
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Model Top  Elevations")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(top,masked_values=masked_values)
cb = plt.colorbar(quadmesh, shrink=0.5, ax=ax)
fig.savefig(os.path.join('fig','top.png'),dpi=300)

# plot bot array
botm = ml.dis.botm.array
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Model Bottom Elevations")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(botm,masked_values=masked_values)
inactive = mapview.plot_inactive()
cb = plt.colorbar(quadmesh, shrink=0.5, ax=ax)
fig.savefig(os.path.join('fig','bot.png'),dpi=300)

# thickness
thk = top - botm
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Model Layer Thickness")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(thk,masked_values=masked_values)
inactive = mapview.plot_inactive()
cb = plt.colorbar(quadmesh, shrink=0.5, ax=ax)
fig.savefig(os.path.join('fig','thk.png'),dpi=300)


'''
# check raster
import rasterio
import gisutils
from gisutils.projection import get_authority_crs
from gisutils.raster import get_values_at_points


(rasterfile, x=None, y=None, band=1,
                         points=None, points_crs=None,
                         xarray_variable=None,
                         out_of_bounds_errors='coerce',
                         method='nearest', size_thresh=1e9):
rasterfile = os.path.join('..','..','gis','dtm_no_drn_ext_trim.tif')
src = rasterio.open(rasterfile)
top = src.read(1)
crs = get_authority_crs(src.crs)

data= get_values_at_points(rasterfile,x,y,points_crs=points_crs)

dis = m.get_package('dis')

x =  m.modelgrid.get_xcellcenters_for_layer(0).ravel()
y =  m.modelgrid.get_ycellcenters_for_layer(0).ravel()

results = src.sample(list(zip(x, y)))
results = np.atleast_1d(np.squeeze(list(results)))
results = results.astype(float)



rasterfile = os.path.join('..','gis','idomain.tif')
idomain = rasterio.open(rasterfile)
top = raster.read(1)
crs = get_authority_crs(raster.crs)

points_crs = get_authority_crs(m.modelgrid.crs)


'''
