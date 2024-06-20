import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import flopy
from mfsetup import MF6model

wd = os.getcwd()

# generate model 
m = MF6model.setup_from_yaml('sescousse.yml',verbose=True)
m.write_input()

# write grid to shapefile 
m.modelgrid.write_shapefile('postproc/shps/grid.shp')

# load simulation and run model 
sim = flopy.mf6.MFSimulation.load(wd)
ml = sim.get_model('sescousse')

# reset some parameters not properly set in yml
sim.ims.outer_dvclose = 1

# rerun simulation
sim.write_simulation()
success, buff = sim.run_simulation()

# plot head 
xlabel = 'x (estward) m RGF93'
ylabel = 'y (northward) m RGF93'
head = flopy.utils.HeadFile('sescousse.hds')
hdata = head.get_alldata()[0]
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Head")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(hdata)#,vmin=19,vmax=21)
cb = plt.colorbar(quadmesh,label='m NGF', shrink=0.5, ax=ax)
fig.tight_layout()

# plot boundary conditions 
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_title("Model Boundary Conditions")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
pmv = flopy.plot.PlotMapView(model=ml, ax=ax)
quadmesh = mapview.plot_ibound()
quadmesh = pmv.plot_bc('riv')
quadmesh = pmv.plot_bc('ghb')
quadmesh = pmv.plot_bc('drn')
fig.tight_layout()

# plot top array
top = ml.dis.top.array
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title("Model Bottom Elevations")
mapview = flopy.plot.PlotMapView(model=ml, layer=0)
quadmesh = mapview.plot_array(top)
inactive = mapview.plot_inactive()
cb = plt.colorbar(quadmesh, shrink=0.5, ax=ax)



