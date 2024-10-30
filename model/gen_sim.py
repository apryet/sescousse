import sys, os, shutil
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy


# load dtm raster
dtm_file = os.path.join('..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)

# template dir 
tpl_dir = 'master_glm'

# --------------------------------------------------------
# no drains  
# ------------------------------------------------------

#  sim dir 
sim_dir = 'nodrn'

# cp simulation from template 
if not os.path.exists(sim_dir):
    shutil.copytree(tpl_dir,sim_dir)
else: 
    print('Failed to copy, directory exists')

# load new model 
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir)
ml = sim.get_model()

# load drain package 
drn = ml.remove_package('drn_0')

# make sure budget and heads are saved
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}

# write new simulation and run 
sim.write_simulation()
success, buff = sim.run_simulation(report=True)

# --------------------------------------------------------
# identical drain network with alternative drainage level 
# ------------------------------------------------------

#  sim dir 
sim_dir = 'drn210'

# drainage depth from surface (dtm) 
dd = 2.1 # 

# cp simulation from template 
if not os.path.exists(sim_dir):
    shutil.copytree(tpl_dir,sim_dir)
else: 
    print('Failed to copy, directory exists')

# load new model 
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir)
ml = sim.get_model()
nper = sim.tdis.nper.data

# resample dtm to model grid 
dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')

# load drain package 
drn = ml.get_package('drn_0')
drn_spd = drn.stress_period_data.get_data()

rec0 = drn_spd[0].copy()

rec0['elev'] = [dtm[x['cellid'][1],x['cellid'][2]]-dd for x in rec0]

# set drain spd 
drn.stress_period_data.set_data({i:rec0 for i in range(nper)}) 

# make sure budget and heads are saved
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}

# write new simulation and run 
sim.write_simulation()
success, buff = sim.run_simulation(report=True)

# ----------------------------------------
# ---- alternative drainage network 
# ----------------------------------------


#  sim dir 
sim_dir = 'drn_110'

# from notebook, deserves sorting 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import shapely
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPoint,
    Point,
    Polygon,
)

import flopy
import flopy.discretization as fgrid
import flopy.plot as fplot
from flopy.utils import GridIntersect


from flopy.utils.geospatial_utils import GeoSpatialCollection


mg = ml.modelgrid
ix = GridIntersect(mg)

drn_lines_shpfile = os.path.join('..','gis','drn_lines.shp')
drn_lines_gsc = GeoSpatialCollection(drn_lines_shp)
drn_lines_geoms = drn_lines_gsc.flopy_geometry
drn_lines_geom = drn_lines_geoms[0]

ls1 = LineString([(95, 105), (30, 50)])
ls2 = LineString([(30, 50), (90, 22)])
ls3 = LineString([(90, 22), (0, 0)])
mls = MultiLineString(lines=[ls1, ls2, ls3])

drn_shp = 
result = ix.intersect(drn_lines_geom)

# plot 
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
mg.plot(ax=ax)
ix.plot_linestring(result, ax=ax, cmap="viridis")

# load initial drn package 

# replace drn package 




