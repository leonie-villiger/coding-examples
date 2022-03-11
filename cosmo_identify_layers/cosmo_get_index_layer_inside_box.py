#!/usr/bin/env python
# coding: utf-8
# ----------------------------------------
# for a series of time steps
#   extract vertical profiles inside lon lat box
#   determine layers cloudtob, cloud base and subcloud
#   calculate cloud fraction at cloud base 
# write time series of layer index, median/mean values and cloud fraction to netcdf
# applicable for output of different COSMOiso simulations
# -----------
# output:
# 1 netcdf (for whole period) with
#   layer index
#   layer median altitude
#   layer median pressure
#   cloudbase cloud fraction
# ----------------------------------------
# Mar 2022 Leonie Villiger
# ----------------------------------------

# ------------------------
# import modules
# ------------------------
import os, sys, glob, gc 
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt

# -------------
# functions and dictionaries
# -------------

# create list of dt.dattime items
def make_datelist(startdate,enddate,hstep):
    nowdate = startdate
    datelist = []
    while nowdate <= enddate:
        datelist.append(nowdate)
        nowdate = nowdate+dt.timedelta(hours=hstep)
    return datelist

# read cosmo iso data inside lon lat box and on specific layer
def read_data_box(ncfile, varname, latidx, lonidx, levidx=None):
    """
    ncfile:  path to ncfile
    varname: name of 4D variable on the ncfile
    latidx:  indices of latitudes inside box
    lonidx:  indices of longitudes inside box
    levidx:  if given; values from this level are read; 
             otherwise the whole profile is taken
    """
    ncfile = Dataset(ncfile)
    var    = np.array(ncfile.variables[varname]) #(time, level, rlat, rlon)
    del ncfile
    gc.collect()
    
    if levidx:
        var = np.squeeze(var[0, levidx, latidx, lonidx])
    else:
        var = np.squeeze(var[0, :, latidx, lonidx])
    return var

# translate var to cosmofile (every simulated time step has 5 output files due to size issues)
cosfiles = {
    'U':'1', 'V':'1', 'W':'1', 'P':'1', 'T':'1', 'POT_VORTIC':'1', #.ncfile1
    'QV':'2', 'QV18O':'2', 'QV2H':'2', #.ncfile2
    'QC':'3', 'QC18O':'3', 'QC2H':'3', 'QI':'3', 'QI18O':'3', 'QI2H':'3', #.ncfile3
    'QR':'4', 'QR18O':'4', 'QR2H':'4', 'QS':'4', 'QS18O':'4', 'QS2H':'4', #.ncfile4
    'PS':'5', 'TOT_PREC':'5', 'TWATER':'5', 'TQV':'5', 'AEVAP_S':'5',
    'TQC':'5', 'TQR':'5' #.ncfile5 -> many more!
    }

# -------------------------------
# setting, ADJUST if needed!

# box to extract data from
latmin,latmax,lonmin,lonmax = 11, 16,-61,-54.5 #precip box
box_name = str(latmin)+'to'+str(latmax)+'N_'+str(lonmin)+'to'+str(lonmax)+'E'

# chose simulation
cos_sim = '1km' #BIG_REF_v3, 5km, 1km
cosname = {'1km': 'COSMO$_\mathrm{iso,1km}$',
           '5km': 'COSMO$_\mathrm{iso,5km}$',
           'BIG_REF_v3': 'COSMO$_\mathrm{iso,10km}$'}
# dirs
idir  = '../cosmo/'+cos_sim+'/'
constfile = idir+'LMCONSTANTS.nc'
ncout = '../cosmo_manipulated/cosmo_'+cos_sim+'_'+box_name+'_ct-cb-sc_alt-cf_new.nc'

# max cloud base height (we are interessted in low-level clouds)
z_top = 1300. #threshold from Vial et al. (2017) Geophys. Rev.
z_top_plt = 3000. # max height of profiles on plot

# thresholds (must exceed) for rain and cloud
th_cloud = 10**(-5) #kgkg-1; threshold from Vial et al. (2017) Geophys. Rev.
th_rain  = 10**(-6) #kgkg-1

# period
start,end = dt.datetime(2020,1,20,0), dt.datetime(2020,2,13,23)
datelist  = make_datelist(startdate=start, enddate=end, hstep=1)
print('\nperiod ', datelist[0], datelist[-1])
ntime     = len(datelist)
print('ntime: ', ntime)

# files
file_list = []
for date in datelist:
    file_list.append(idir+'lffd'+date.strftime('%Y%m%d%H'))

# varlist 
varlist = ['P', 'QC']  

# -------------------------------
# get grid from cosmo
coords    = Dataset(constfile,'r')
hhl       = coords.variables['HHL'][0,:,:,:]   #altitude -> middle of layer
Z3D       = (hhl[1:, ...] + hhl[:-1, ...]) / 2 #compute z coordinates from altitude
z_n       = Z3D.shape[0]
lon       = coords.variables['lon'][:]
lat       = coords.variables['lat'][:]
del coords
gc.collect()
print('\nn levels: ', z_n)
print('lowest model level: ', Z3D[0,:,:])
print('highest model level: ', Z3D[-1,:,:])

# --------------------------------
# find indixes of grid cell inside box
latidx,lonidx = np.where((lon>=lonmin)&(lon<=lonmax)&(lat>=latmin)&(lat<=latmax))
npos    = len(latidx)
print('\nnumber of grid cells inside box: ', npos)

# -------------------------------
# create empty objects to store output of loop
timecoord      = np.zeros(0) # time dimension, datetime format
timecoord_nc   = np.zeros(0) # time dimension, sec since format
out = dict()
na_array = np.empty((ntime))
na_array[:] = np.nan
out['cloud_n']    = na_array.copy() #number of positions with liquid clouds (used for mean calc)
out['cb_idx']     = na_array.copy() #median index for layer
out['sc_idx']     = na_array.copy() 
out['ct_idx']     = na_array.copy()
out['med_Z_cb']   = na_array.copy() #median altitude for layer
out['med_Z_sc']   = na_array.copy()
out['med_Z_ct']   = na_array.copy()
out['med_P_cb']   = na_array.copy() #median pressure for layer
out['med_P_sc']   = na_array.copy()
out['med_P_ct']   = na_array.copy()
del na_array

# -------------------------------
# extract profiles and find atmospheric layers 
# store index, altitude, pressure in netcdf

# run through files
for ts in range(ntime):
    
    # -------------------------------------------------
    # read data 
    
    # read file
    ifile = file_list[ts]
    print('\nfile ', ifile)

    # get time step
    date_dt      = dt.datetime(int(ifile[-10:-6]),int(ifile[-6:-4]),int(ifile[-4:-2]),int(ifile[-2:]))
    timecoord    = np.append(timecoord,date_dt)
    timecoord_nc = np.append(timecoord_nc, (date_dt-dt.datetime(2020,1,1,0,0,0)).total_seconds())
    print('time step: ', date_dt)
    del date_dt

    # object to store temporary data
    dbox = {}

    # altitude profile inside box
    dbox['Z'] = np.squeeze(Z3D[:, latidx, lonidx]).T #(level, rlat, rlon)
    print('\nshape of Z in box ', dbox['Z'].shape)

    # read data from file
    for var in varlist:
        dbox[var] = read_data_box(ncfile=ifile+'.ncfile'+cosfiles[var], varname=var,latidx=latidx, lonidx=lonidx, levidx=None)
        print('shape of ',var,' in box: ', dbox[var].shape)
    del ifile

    # find profiles with low level liquid cloud --> used to identify layers
    print('\nfind layers')

    # prepare objects to store identified indices
    idxs_cloud  = [] #profiles with low level clouds
    idxs_cb     = [] #cloud base
    idxs_sc     = [] #subcloud layer
    idxs_ct     = [] #cloud top

    # assign profiles to categories
    for pos in range(npos):

        # get profiles
        Z  = np.squeeze(dbox['Z'][pos,:]) #whole profile
        QC = np.squeeze(dbox['QC'][pos,:]) #kg/kg 

        # check if low level liquid clouds are present
        if np.sum(QC[Z<z_top]>th_cloud)>0: 
            
            # profiles with liquid clouds
            idxs_cloud.append(pos)
             
            # find CLOUD BASE - first level from surface upwards with liquid cloud
            idxs_cb.append(np.where(QC>th_cloud)[0][-1])
            
            # find SUBCLOUD - mid index between surface and cloud base
            idxs_sc.append(int(np.round((idxs_cb[-1] + (z_n-1))*0.5))) 
            
            # find CLOUD TOP - start at grid cell above cloud base and go upwards until first grid cell without liquid cloud is found (note vertical dimension in COSMO is reversed)
            for k in np.arange(idxs_cb[-1]-1,0,-1):
                if QC[k]==0:
                    idxs_ct.append(int(k))
                    break
    
    # remove objects no longer needed 
    del QC,Z

    # -------------------------------------------------
    # if any clouds: identify median indices for layers 
    print('\nstore identified layers')
    
    if len(idxs_cloud)>0:
        # get number of clouds and indices
        out['cloud_n'][ts] = len(idxs_cloud) 
        out['cb_idx'][ts]  = np.round(np.median(idxs_cb))
        out['sc_idx'][ts]  = np.round(np.median(idxs_sc))
        out['ct_idx'][ts]  = np.round(np.median(idxs_ct))

        # get median altitude/pressure
        out['med_Z_cb'][ts] = np.median(dbox['Z'][:,int(out['cb_idx'][ts])])
        out['med_P_cb'][ts] = np.median(dbox['P'][:,int(out['cb_idx'][ts])])
        out['med_Z_sc'][ts] = np.median(dbox['Z'][:,int(out['sc_idx'][ts])])        
        out['med_P_sc'][ts] = np.median(dbox['P'][:,int(out['sc_idx'][ts])])
        out['med_Z_ct'][ts] = np.median(dbox['Z'][:,int(out['ct_idx'][ts])])
        out['med_P_ct'][ts] = np.median(dbox['P'][:,int(out['ct_idx'][ts])])

    else:
        out['cloud_n'][ts] = 0
        print('No clouds in box. Returning NaN values.')

    print('profiles with clouds:  ', out['cloud_n'][ts])
    print('cloud top, height: ', out['med_Z_ct'][ts])
    print('cloud base, height: ', out['med_Z_cb'][ts])
    print('subcloud, height:  ', out['med_Z_sc'][ts])
    
    # remove objects/clear memory
    del dbox, idxs_cloud, idxs_cb, idxs_sc, idxs_ct
    gc.collect()

# ----------------------
print('')
print('store netcdf')

# --- open a netCDF file to write
with Dataset(ncout, mode='w') as ncfile:

    # --- global attributes
    ncfile.description = 'COSMOiso ('+cos_sim+') index, median of Z and P for cloud top, cloud base and subcloud layer, and cloud fraction at cloud base (grid points with QC>'+str(th_cloud)+'kgkg-1) at cloudbase inside lon-lat box ('+str(lonmin)+'°E, '+str(latmin)+'°N to '+str(lonmax)+'°E, '+str(latmax)+'°N'+'), time period: '+start.strftime('%Y%m%d_%H')+'-'+end.strftime('%Y%m%d_%H')
    ncfile.history = 'Created in Zürich, '+dt.date.today().strftime('%Y-%m-%d')
    
    # --- create dimensions
    ncfile.createDimension('time', ntime)
    ncfile.createDimension('grid_points_inside_box', npos)

    # --- write variables
    time                = ncfile.createVariable('time', 'f8', ('time'))
    time.units          = 'seconds since 2020-01-01 00:00:00'
    time.long_name      = 'time'
    ncfile['time'][:]   = timecoord_nc

    liqcld_n               = ncfile.createVariable('liqcld_n', 'f8', ('time'))
    liqcld_n.units         = 'count'
    liqcld_n.long_name     = 'number of positions with a liquid cloud (QV>'+str(th_cloud)+'kgkg-1) detected.'
    ncfile['liqcld_n'][:]  = out['cloud_n']

    cf               = ncfile.createVariable('cf', 'f8', ('time'))
    cf.units         = 'ratio'
    cf.long_name     = 'cloud fraction at cloud base (grid points with QV>'+str(th_cloud)+'kgkg-1)'
    ncfile['cf'][:]  = out['cloud_n']/npos

    idx_ct              = ncfile.createVariable('idx_ct', 'f8', ('time'))
    idx_ct.units        = 'index'
    idx_ct.long_name    = 'median index of cloud top'
    ncfile['idx_ct'][:] = out['ct_idx']

    idx_cb              = ncfile.createVariable('idx_cb', 'f8', ('time'))
    idx_cb.units        = 'index'
    idx_cb.long_name    = 'median index of cloud base'
    ncfile['idx_cb'][:] = out['cb_idx']

    idx_sc              = ncfile.createVariable('idx_sc', 'f8', ('time'))
    idx_sc.units        = 'index'
    idx_sc.long_name    = 'median index of subcloud'
    ncfile['idx_sc'][:] = out['sc_idx']

    Z_ct              = ncfile.createVariable('Z_ct', 'f8', ('time'))
    Z_ct.units        = 'm'
    Z_ct.long_name    = 'median altitude of cloud top'
    ncfile['Z_ct'][:] = out['med_Z_ct']

    Z_cb              = ncfile.createVariable('Z_cb', 'f8', ('time'))
    Z_cb.units        = 'm'
    Z_cb.long_name    = 'median altitude of cloud base'
    ncfile['Z_cb'][:] = out['med_Z_cb']

    Z_sc              = ncfile.createVariable('Z_sc', 'f8', ('time'))
    Z_sc.units        = 'm'
    Z_sc.long_name    = 'median altitude of subcloud'
    ncfile['Z_sc'][:] = out['med_Z_sc']

    P_ct              = ncfile.createVariable('P_ct', 'f8', ('time'))
    P_ct.units        = 'Pa'
    P_ct.long_name    = 'median pressure of cloud top'
    ncfile['P_ct'][:] = out['med_P_ct']

    P_cb              = ncfile.createVariable('P_cb', 'f8', ('time'))
    P_cb.units        = 'Pa'
    P_cb.long_name    = 'median pressure of cloud base'
    ncfile['P_cb'][:] = out['med_P_cb']

    P_sc              = ncfile.createVariable('P_sc', 'f8', ('time'))
    P_sc.units        = 'Pa'
    P_sc.long_name    = 'median pressure of subcloud'
    ncfile['P_sc'][:] = out['med_P_sc']

