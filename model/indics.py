import os
import numpy as np
import pandas as pd
import flopy

print('Loading rasters')
# load idomain raster
idomain_file = os.path.join('..','gis','idomain.tif') 
idomain_rast = flopy.utils.Raster.load(idomain_file)
# load dtm raster
dtm_file = os.path.join('..','gis','dtm_no_drn_ext_trim_filt.tif') 
dtm_rast = flopy.utils.Raster.load(dtm_file)

# -------------------------------------
# ----  indicator definitions     -----
# -------------------------------------

def get_indic(sim_dir):
    # alert thresholds (depths)
    depth_w = 0.4 # water excess 
    depth_d = 1.2 # water deficit
    # load mf model and simulated heads
    # skip all packages but tdis to save time (otherwise takes ages!)
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_dir,load_only=['tdis'],verbosity_level=0)
    ml = sim.get_model()
    nper = sim.tdis.nper.data
    headfile = flopy.utils.binaryfile.HeadFile(os.path.join(sim_dir,'sescousse.hds'))
    hds = headfile.get_alldata()
    # resample idomain and dtm 
    idomain = idomain_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
    idomain_3d = np.stack([idomain]*nper) # transient idomain 
    dtm = dtm_rast.resample_to_grid(ml.modelgrid, band=1, method='nearest')
    # number of cells in domain of interest (idomain==3)
    ncells = (idomain==3).sum()
    # critical levels 
    z_w = dtm - depth_w # m NGF
    z_d = dtm - depth_d # m NGF
    # dates out 
    start_date= pd.to_datetime(sim.tdis.start_date_time.get_data())
    end_date = pd.to_datetime(start_date+ pd.to_timedelta(nper-1,'d'))
    dates_out  = pd.date_range(start_date,end_date).date
    # masked head array out of area of interest (idomain==3)
    mhds = np.ma.masked_where(idomain_3d<3, hds[:,0,:,:])
    # records of spatially averaged water excess/stress
    w_records = ((mhds-z_w)*(mhds>z_w)).sum(axis=(1,2))/ncells
    d_records = ((mhds-z_d)*(mhds<z_d)).sum(axis=(1,2))/ncells
    # aggregate in df
    df = pd.DataFrame({'w':w_records,'d':d_records},index=dates_out)
    return(df)


if __name__ == "__main__":
    print('Loading rasters')
    # load idomain raster
    idomain_file = os.path.join('..','..','gis','idomain.tif') 
    idomain_rast = flopy.utils.Raster.load(idomain_file)
    # load dtm raster
    dtm_file = os.path.join('..','..','gis','dtm_no_drn_ext_trim_filt.tif') 
    dtm_rast = flopy.utils.Raster.load(dtm_file)
    # compute indic 
    print('Processing indicators')
    df = get_indic('.')
    # write csv
    df.to_csv('indics.csv')
    print('Indicators saved to indics.csv')

