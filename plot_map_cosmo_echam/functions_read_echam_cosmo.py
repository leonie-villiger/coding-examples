#!/usr/bin/env python
# coding: utf-8
# ----------------------------------------
# read data from ECHAMwiso, COSMO 
# interpolate to certain pressure level
# ----------------------------------------
# Jan 2022 Leonie Villiger
# ----------------------------------------

# ------------------------
# import modules
# ------------------------
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import functions_netcdf as fn
import metpy.calc
from metpy.units import units

# ------------------------
# read ECHAMwiso
# ------------------------
def read_echam(idir_echam, var, timestep, lev):
    """
    idir_echam: path to cfsd* files
    var:        variable to read/calc (only a bunch impelemnted in function)
    timestep:   time step to read
    lev:        pressure level in hPa to interpolate to
    """
    # prepare output
    echam_data = dict()

    # read file
    print('reading file: ', idir_echam+'cfsd'+timestep.strftime('%Y%m%d%H')+'.nc')
    ncfile = Dataset(idir_echam+'cfsd'+timestep.strftime('%Y%m%d%H')+'.nc')

    # horizontal grid
    lon = np.array(ncfile.variables['lon'])
    lat = np.array(ncfile.variables['lat'])
    longrid, latgrid = np.meshgrid(lon,lat)
    echam_data['lon'] = longrid.copy()
    echam_data['lat'] = latgrid.copy()
    del lon, lat, longrid, latgrid

    # prepare interpolation
    #nf.printidim(ncfile)
    ps    = np.squeeze(np.array(ncfile.variables['PS'][0,:,:]))/100.
    nlev  = fn.readcdfdim(idir_echam+'cfsd'+timestep.strftime('%Y%m%d%H')+'.nc','level')
    ak,bk = fn.echamakbk(idir_echam+'cfsd'+timestep.strftime('%Y%m%d%H')+'.nc')
    p3d   = fn.getP3D(ps,ak,bk,nlev)
    del nlev, ps, ak, bk

    # read variables and interpolate to pressure level
    if var in ['QV', 'QC']:
        var3d = np.squeeze(np.array(ncfile.variables[var][0,:,:,:]))*1000.
        echam_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del var3d

    elif var=='dD':
        qv2h  = np.squeeze(np.array(ncfile.variables['QV2H'][0,:,:,:]))
        var3d = (qv2h-1)*1000.
        echam_data['dD'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv2h, var3d

    elif var=='d18O':
        qv18o = np.squeeze(np.array(ncfile.variables['QV18O'][0,:,:,:]))
        var3d = (qv18o-1)*1000.
        echam_data['d18O'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv18o,var3d

    elif var=='dexc':
        qv2h  = np.squeeze(np.array(ncfile.variables['QV2H'][0,:,:,:]))
        qv18o = np.squeeze(np.array(ncfile.variables['QV18O'][0,:,:,:]))
        var3d = ((qv2h-1)*1000.)-8*((qv18o-1)*1000.)
        echam_data['dexc'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv2h, qv18o, var3d
    
    elif var=='T': #degK->degC
        var3d = np.squeeze(np.array(ncfile.variables[var][0,:,:,:]))-273.15
        echam_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del var3d

    elif var=='virtpotT':
        T     = np.squeeze(np.array(ncfile.variables['T'][0,:,:,:])) #Kelvin
        q     = np.squeeze(np.array(ncfile.variables['QV'][0,:,:,:])) #kgkg-1
        mr    = metpy.calc.mixing_ratio_from_specific_humidity(specific_humidity=q*units('kg/kg'))
        mr    = np.array(mr) #kgkg-1
        var3d = metpy.calc.virtual_potential_temperature(pressure=p3d*units('hPa'), 
            temperature=T*units('kelvin'), mixing=mr*units('kg/kg'))
        var3d = np.array(var3d)-273.15 #kelvin->degC
        echam_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del T, q, mr, var3d
    else:
        print('implement code to read '+var+' from ECHAMwiso file')

    # delete stuff no longer needed
    del ncfile, p3d

    return echam_data

# ------------------------
# read COSMO data
# ------------------------
def read_cosmo(idir_cosmo, var, timestep, lev):
    """
    idir_echam: path to cfsd* files
    var:        variable to read/calc (only a bunch impelemnted in function)
    timestep:   time step to read
    lev:        pressure level in hPa to interpolate to
    """

    # translate var to cosmofile
    cosfiles = {
        'U':'1', 'V':'1', 'W':'1', 'P':'1', 'T':'1', 'POT_VORTIC':'1', #.ncfile1
        'QV':'2', 'QV18O':'2', 'QV2H':'2', #.ncfile2
        'QC':'3', 'QC18O':'3', 'QC2H':'3', 'QI':'3', 'QI18O':'3', 'QI2H':'3', #.ncfile3
        'QR':'4', 'QR18O':'4', 'QR2H':'4', 'QS':'4', 'QS18O':'4', 'QS2H':'4', #.ncfile4
        'PS':'5', 'TOT_PREC':'5', 'TWATER':'5', 'TQV':'5', 'AEVAP_S':'5',
        'TQC':'5', 'TQR':'5' #.ncfile5 -> many more!
        }

    # prepare output
    cosmo_data = dict()

    # horizontal grid 
    coords  = Dataset(idir_cosmo+'LMCONSTANTS.nc','r')
    longrid = coords.variables['lon'][:]
    latgrid = coords.variables['lat'][:]
    cosmo_data['lon'] = longrid.copy()
    cosmo_data['lat'] = latgrid.copy()
    del coords, longrid, latgrid

    # prepare interpolation
    print('reading files: ',idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile*')
    ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile1')
    p3d    = np.squeeze(np.array(ncfile.variables['P'][0,:,:,:]))/100. #Pa->hPa
    del ncfile

    # read variable and interpolate to pressure level
    if var in ['QV', 'QR', 'QC']: #variables scales by 1000. kg/kg -> g/kg
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles[var])
        var3d  = np.squeeze(np.array(ncfile.variables[var][0,:,:,:]))*1000.
        cosmo_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del var3d, p3d, ncfile

    elif var=='dD':
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles['QV'])
        qv     = np.squeeze(np.array(ncfile.variables['QV'][0,:,:,:]))
        qv2h   = np.squeeze(np.array(ncfile.variables['QV2H'][0,:,:,:]))
        var3d  = (qv2h/qv-1)*1000.
        cosmo_data['dD'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv, qv2h, var3d, p3d, ncfile

    elif var=='d18O':
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles['QV'])
        qv     = np.squeeze(np.array(ncfile.variables['QV'][0,:,:,:]))
        qv18o  = np.squeeze(np.array(ncfile.variables['QV18O'][0,:,:,:]))
        var3d  = (qv18o/qv-1)*1000.
        cosmo_data['d18O'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv, qv18o, var3d, p3d, ncfile

    elif var=='dexc':
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles['QV'])
        qv     = np.squeeze(np.array(ncfile.variables['QV'][0,:,:,:]))
        qv2h   = np.squeeze(np.array(ncfile.variables['QV2H'][0,:,:,:]))
        qv18o  = np.squeeze(np.array(ncfile.variables['QV18O'][0,:,:,:]))
        var3d  = ((qv2h/qv-1)*1000.)-8*((qv18o/qv-1)*1000.)
        cosmo_data['dexc'] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del qv, qv2h, qv18o, var3d, p3d, ncfile

    elif var=='T': #degK->degC
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles[var])
        var3d  = np.squeeze(np.array(ncfile.variables[var][0,:,:,:]))-273.15
        cosmo_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del var3d, p3d, ncfile

    elif var=='virtpotT':
        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles['T'])
        T     = np.squeeze(np.array(ncfile.variables['T'][0,:,:,:])) #Kelvin
        del ncfile

        ncfile = Dataset(idir_cosmo+'lffd'+timestep.strftime('%Y%m%d%H')+'.ncfile'+cosfiles['QV'])
        q     = np.squeeze(np.array(ncfile.variables['QV'][0,:,:,:])) #kgkg-1
        del ncfile

        mr    = metpy.calc.mixing_ratio_from_specific_humidity(specific_humidity=q*units('kg/kg'))
        mr    = np.array(mr) #kgkg-1
        
        var3d = metpy.calc.virtual_potential_temperature(pressure=p3d*units('hPa'),          
            temperature=T*units('kelvin'), mixing=mr*units('kg/kg'))
        var3d = np.array(var3d)-273.15 #kelvin->degC
        cosmo_data[var] = fn.inter2lev(var3d, p3d, int(lev), mode='p')
        del T, q, mr, var3d

    else:
        print('implement code to read '+var+' from COSMOiso file')

    return cosmo_data
