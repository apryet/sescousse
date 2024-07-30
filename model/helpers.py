
def process_secondary_obs(ws='.'):
    # load dependencies inside the function so that they get carried over to forward_run.py by PstFrom
    import os
    import pandas as pd

    def write_tdif_obs(orgf, newf, ws='.'):
        df = pd.read_csv(os.path.join(ws,orgf), index_col='time')
        df = df - df.mean()
        df.to_csv(os.path.join(ws,newf))
        return

    # write the tdiff observation csv's
    write_tdif_obs('sescousse.head.csv', 'sescousse.head.tdiff.csv', ws)

    print('Secondary observation files processed.')
    return

def set_mf_par_vals():
    ''' 
    Set Modflow parameter values 
    '''
    # load dependencies inside the function so that they get carried over to forward_run.py by PstFrom
    import os
    import pandas as pd
    import flopy

    # load model with flopy
    print('Loading model...')
    parvals = pd.read_csv('par.dat',header=None,delim_whitespace=True,index_col=0)[1]
    sim = flopy.mf6.MFSimulation.load(sim_ws='.')
    ml = sim.get_model()

    # load update swb vars 
    # recharge (soil drainage) and aquifer root water uptake
    swb = pd.read_csv('swb_vars.csv',index_col=0,parse_dates=True)
    rech = swb.loc[:,'R']*0.001/86400 # mm/d to m/s
    evt = swb.loc[:,'RU']*0.001/86400 # mm/d to m/s

    # set initial values for steady state
    rech.iloc[0] = 0. # 9.5e-9  # 300mm/year in m/s, for initial steady state
    evt.iloc[0] = 0. # deactivate evt for initial steady state

    # set recharge 
    print('Setting recharge...')
    rcha = ml.get_package('rcha_0')
    rcha.recharge = {i:rech[i] for i in range(ml.nper)}

    # set transpiration
    print('Setting evt...')
    evta = ml.get_package('evt')
    evta.rate = {i:evt[i] for i in range(ml.nper)}

    # set k 
    print('Setting k values...')
    ml.npf.k.set_data(parvals.loc['k'])

    # set sy 
    print('Setting sy values...')
    ml.sto.sy.set_data(parvals.loc['sy'])

    # set cghb
    print('Setting cghb values...')
    ghb_spd = ml.ghb.stress_period_data.get_data()
    for i in ghb_spd.keys():
        ghb_rec = ghb_spd[i]
        # western bc
        idx = ghb_rec['boundname']=='w'
        ghb_rec['cond'][idx] = parvals.loc['cghbw']
        # eastern bc
        idx = ghb_rec['boundname']=='e'
        ghb_rec['cond'][idx] = parvals.loc['cghbe']

    # update values 
    ml.ghb.stress_period_data.set_data(ghb_spd)

    # set cdrn
    print('Setting cdrn values...')
    drn = ml.get_package('drn_0')
    drn_spd = drn.stress_period_data.get_data()
    drn_spd[0]['cond']=parvals.loc['cdrn']
    drn.stress_period_data.set_data(drn_spd)

    # set criv
    print('Setting criv values...')
    riv = ml.get_package('riv')
    riv_spd = riv.stress_period_data.get_data()
    for i in riv_spd.keys():
        riv_spd[i]['cond'] = parvals.loc['criv']

    riv.stress_period_data.set_data(riv_spd)

    print('Writing new model files...')
    sim.write_simulation()


def write_tot_drn(ws='.'):
    import os 
    import pandas as pd 
    zone_surf = 1071422. # model domain surface (area of interest only), m2
    drnobs = pd.read_csv(os.path.join(ws,'sescousse.drn.obs.output.csv'),index_col=0)
    drnobs['flow']= drnobs.DRN*(-1/zone_surf*1000*86400) # m3/s to mm/d
    tot_drn = drnobs.flow.sum()
    tot_drn_df = pd.DataFrame({'tot_drn_mm':[tot_drn]})
    tot_drn_df.index.name = 'time'
    tot_drn_df.to_csv(os.path.join(ws,'tot_drn.csv'))


