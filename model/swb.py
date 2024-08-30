import os
import numpy as np
import pandas as pd

# modified from :
# https://github.com/soilwater/pynotes-agriscience/blob/gh-pages/exercises/soil_water_balance.ipynb
# removed runoff, set Kc to constant, external computation of ET0
def run_swb(clim_file='clim.csv',theta_sat = None, D_max= None, par_file=None, cwd='.'):
    # load clim data
    clim = pd.read_csv(os.path.join(cwd,clim_file),index_col=0, parse_dates=True)
    P = clim.P
    PET = clim.PET

    # set start/end dates
    start_date = clim.index.min()
    end_date   = clim.index.max()
    sim_dates = pd.date_range(start_date,end_date) 

    # vegetation and soil parameters
    #I_max = 1 # Canopy and crop residue interception (evaporation) in mm/day
    z = 400 # Length of the soil profile in mm

    # if par_file provided, read file and set par values 
    if par_file is not None:
        # load par values 
        parvals = pd.read_csv('par.dat',header=None,delim_whitespace=True,index_col=0)[1]
        # set parameter values 
        theta_sat = parvals.loc['tsat']
        D_max = parvals.loc['dmax']

    # characteristic water contents
    if theta_sat is None: 
        theta_sat = 0.40

    # Maximum daily drainage rate [mm/d]
    if D_max is None: 
        D_max = 25 

    # with adjustable theta_sat, fc and r have to be defined as ratios to avoid conflicts
    # (estimates deserve to be refined)
    theta_fc = 0.5*theta_sat
    theta_r  = 0.2*theta_sat

    # crop coefficient 
    Kc = 1.0 # REF? 

    S_max =theta_sat*z # Saturation
    FC = theta_fc*z  # Field capacity
    PWP = theta_r*z # Permanent Wilting Point

    # allocate arrays with NaN values
    N = sim_dates.shape[0]
    T = np.ones(N)*np.nan # plant transpiration
    D = np.ones(N)*np.nan # soil drainage
    F = np.ones(N)*np.nan # soil fast drainage (overflow)
    S = np.ones(N)*np.nan # soil water storage 
    R = np.ones(N)*np.nan # groundwater recharge

    # Define initial conditions
    S[0] = FC # field capacity
    T[0] = PET[0]*Kc # transpiration from soil
    D[0] = D_max*(S[0]/S_max)**((D_max/S_max - 1)/(-D_max/S_max)) # soil drainage
    F[0] = np.maximum(S[0]-S_max,0) # soil fast drainage
    R[0] = D[0] + F[0] # groundwater recharge

    # Model
    for t in range(1, N):
        # drainage (derived from Wilcox ?)
        D[t] = D_max*(S[t-1]/S_max)**((D_max/S_max - 1)/(-D_max/S_max)) 
        # plant available water
        PAW = np.maximum((S[t-1] - PWP) / (FC - PWP), 0)
        # stress coefficient
        Ks = np.maximum(1 - np.exp(-10*PAW), 0)
        # transpiration 
        T[t] = PET[t] * Kc * Ks
        # update soil water content 
        S[t] = max(S[t-1] + P[t] - T[t] - D[t],0)
        D[t] = S[t-1] - S[t] + P[t] - T[t] 
        # overflow (added to drainage) 
        F[t] = np.maximum(S[t]-S_max,0)
        S[t] = S[t] - F[t]
        R[t] = D[t] + F[t]

    # save recharge and transpiration for modflow  
    swb_df = pd.DataFrame({'P':P, 'PET':PET, 'T':T, 'D':D,'R':R, 'F':F,'S':S},index=sim_dates)

    # consider aquifer root water uptake as evapotranspiration deficit 
    # (water that cannot be taken from the soil is taken from the aquifer)
    swb_df['RU'] = swb_df['PET'] - swb_df['T']

    swb_df.index.name = 'date'
    swb_df.to_csv(os.path.join(cwd,'swb_vars.csv'))

    # save cumulative flows (to be used as predictions)
    swb_cum = pd.DataFrame(swb_df[['T','R','RU']].sum()).T
    swb_cum.index.name = 'time'
    swb_cum.to_csv(os.path.join(cwd,'swb_cum.csv'))



