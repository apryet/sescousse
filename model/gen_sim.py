import os, shutil
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import flopy


# template dir 
tpl_dir = 'master_glm'



# alternative model with deeper drainage level 

# change in drainage level 
dh = -(1.10 - 0.75) # -(new depth - former depth)

#  sim dir 
sim_dir = 'drn110'

# cp 
if not os.path.exists(sim_dir):
    shutil.copytree(tpl_dir,sim_dir)
else: 
    print('Failed to copy, directory exists')


# load new model 
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir)
ml = sim.get_model()
nper = sim.tdis.nper.data

drn = ml.get_package('drn_0')
drn_spd = drn.stress_period_data.get_data()

new_drn_spd = {}

for i,drn_rec in drn_spd.items():
    new_drn_rec = drn_rec.copy()
    new_drn_rec['elev'] = drn_rec['elev'] + dh
    new_drn_spd[i]=new_drn_rec

# set drainage level 
drn_spd = drn.stress_period_data.set_data(new_drn_spd)

# make sure budget and heads are saved
ml.oc.saverecord = {k:[('HEAD','LAST'), ('BUDGET','LAST')]for k in range(nper)}

# write new simulation and run 
sim.write_simulation()
success, buff = sim.run_simulation(report=True)



