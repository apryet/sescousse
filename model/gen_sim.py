import sys, os, shutil, glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy
import shapefile

# load dtm raster
dtm_file = os.path.join('..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)

# template dir 
tpl_dir = 'master_glm'

# cp simulation from template 
def cp_from_tpl(tpl_dir,sim_dir):
    if not os.path.exists(sim_dir):
        shutil.copytree(tpl_dir,sim_dir)
        # clear calibration and output files
        for pattern in ['cal*','*.jcb','*.ins','*.tpl',\
                '*ut.csv','*.list','*.hds','*.cbc']:
            for f in glob.glob(os.path.join(sim_dir,pattern)):
                os.remove(f)
        # clear figures
        for f in glob.glob(os.path.join(sim_dir,'fig','*.*')):
            os.remove(f)
    else: 
        print('WARNING : Failed to copy, directory exists')

# --------------------------------------------------------
# no drains  
# ------------------------------------------------------

#  sim dir 
sim_dir = 'drn0'

# cp simulation from template 
cp_from_tpl(tpl_dir,sim_dir)

# load template model 
sim = flopy.mf6.MFSimulation.load(sim_ws=tpl_dir)
nper = sim.tdis.nper.data
ml = sim.get_model()

# load drain package 
drn = ml.remove_package('drn_0')

# make sure budget and heads are saved
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}

# write new simulation and run
sim.set_sim_path(sim_dir)
sim.write_simulation()
#success, buff = sim.run_simulation(report=True)

# --------------------------------------------------------
# identical drain network with alternative drainage level 
# ------------------------------------------------------

# generate simulation with drainage depth dd (m from surface)
def set_dd(dd,tpl_dir):

    #  sim dir 
    sim_dir = f'drn{int(dd*100)}'

    # cp simulation from template 
    cp_from_tpl(tpl_dir,sim_dir)

    # load new model 
    sim = flopy.mf6.MFSimulation.load(sim_ws=tpl_dir)
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
    ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')] for k in range(nper)}

    # write new simulation and run 
    sim.set_sim_path(sim_dir)
    sim.write_simulation()
    #success, buff = sim.run_simulation(report=True)



for dd in [0.4,1.10,1.50,2.0,3.0]: set_dd(dd,tpl_dir)


# ----------------------------------------
# ---- alternative drainage network 
# ----------------------------------------

# generate simulation with alternative drainage network (dn) and depth (dd, m from land surface)
def set_dn(drn_shp_name, dd, tpl_dir):
    '''
    drn_shp_name : name of shapefile, to be found in gis dir
    dd : drainage depth
    '''
    #  sim dir 
    dn_suffix = drn_shp_name.split('.')[0].split('_')[-1] # yields "hd" or "ld"
    sim_dir = f'drn_alt_{dn_suffix}_{int(dd*100)}' # yields e.g. "drn_alt_ld_110"
    # cp simulation from template 
    cp_from_tpl(tpl_dir,sim_dir)
    # load new model 
    sim = flopy.mf6.MFSimulation.load(sim_ws=tpl_dir)
    ml = sim.get_model()
    nper = sim.tdis.nper.data
    # resample dtm to model grid 
    dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
    # load alternative drainage network
    mg = ml.modelgrid
    ix = flopy.utils.GridIntersect(mg,method='vertex')
    drn_lines_shp = os.path.join('..','gis',drn_shp_name)
    sf = shapefile.Reader(drn_lines_shp)
    shapes = sf.shapes()
    cellids = [x['cellids'] for s in shapes for x in ix.intersects(s)]
    rows = np.array([int(c[0]) for c in cellids]).astype(int)
    cols = np.array([int(c[1]) for c in cellids]).astype(int)
    # load original drain package 
    drn = ml.get_package('drn_0')
    drn_spd = drn.stress_period_data.get_data()
    rec0 = drn_spd[0].copy()
    cdrn = float(rec0['cond'][0])
    # initialize new spd reccarray
    rec0 = flopy.mf6.ModflowGwfdrn.stress_period_data.empty(ml,maxbound=len(rows),boundnames=True)[0]
    rec0['cellid'] = [ (0,int(row),int(col)) for row,col in zip(rows,cols)]
    rec0['elev'] = [dtm[x['cellid'][1],x['cellid'][2]]-dd for x in rec0]
    rec0['cond'] = cdrn
    rec0['boundname'] = 'D' # arbitrary 
    # reset drainage spd
    drn.stress_period_data.set_data({i:rec0 for i in range(nper)}) 
    drn.maxbound = len(rows)
    # update obs package 
    # build obs data
    drn_obs = {'sescousse.drn.obs.output.csv': [('D', 'DRN', 'D')]}
    # initialize obs package
    drn.obs.initialize(filename='sescousse.drn.obs', continuous=drn_obs)
    # make sure budget and heads are saved
    ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')] for k in range(nper)}
    # write new simulation and run
    sim.set_sim_path(sim_dir)
    sim.write_simulation()
    #success, buff = sim.run_simulation(report=True)


# drainage depths
dds = [0.4,1.10]*2

# drainage networks in gis dir
drn_shp_names = ['drn_lines_hd.shp']*2+['drn_lines_ld.shp']*2

for drn_shp_name,dd in zip(drn_shp_names,dds): set_dn(drn_shp_name, dd, tpl_dir)















