#!/usr/bin/python
#-*- coding:utf-8 -*-
# ---------------------------------------
# functions for netCDF handling in python
# ---------------------------------------

# modules
import numpy as np
from netCDF4 import Dataset
from dypy.small_tools import interpolate

# functions
def printivar(ncfile):
    '''
    print info of variables in ncfile
    '''
    infile = Dataset(ncfile, mode='r')
    nc_vars = [var for var in infile.variables]  # list of nc variables
    print('NetCDF variable information:')
    for var in nc_vars:
        print('\tName: ', var)
        print('\t\tdimensions: ', infile.variables[var].dimensions)
        print('\t\tsize: ', infile.variables[var].size)
        
def printidim(ncfile):
    '''
    print info of dimensions in ncfile
    '''
    infile = Dataset(ncfile, mode='r')
    nc_dims = [dim for dim in infile.dimensions]  # list of nc variables
    print('NetCDF dimension information:')
    for dim in nc_dims:
        print('\tName: ', dim)
        print('\t\tsize: ', infile.dimensions[dim].size)

def readcdf(ncfile,varnam):
    '''
    read variable from ncfile
    '''
    infile = Dataset(ncfile, mode='r')
    var = np.array(infile.variables[varnam][:])
    return(var)

def readcdfdim(ncfile,dimnam):
    '''
    read dimension size from ncfile
    '''
    infile = Dataset(ncfile, mode='r')
    dim = int(infile.dimensions[dimnam].size)
    return(dim)

def getP3D(ps,ak,bk,nlev):
    '''
    get pressure for each model level
    ---
    ps: surface pressure (2D: lon,lat)
    ak,bk: aklay, bklay from constant file (1D: nlev)
    nlev: number of vertical levels (integer)
    '''
    ps3d = np.tile(ps,(nlev,1,1))
    p3d = (ak/100.+bk*ps3d.T).T
    return(p3d)

def inter2lev(var3D, lev3D, lev, mode):
    '''
    Interpolates 3-D (level, lat, lon) over level for variable array varr with
    associated pressure grid parr to the scalar pressure level plevel
    ---
    var3D: var to interpolate (lev, lat, lon)
    lev3D: level information, e.g. pressure, theta
    lev: level to interpolate to, e.g. pressure, theta value
    mode: p for pressure, th for theta
    NOTE: you need to check with dataset, whether order of var3D, lev3D is correct
    '''
    if mode=='p':
        v_i = interpolate(var3D[::1,:,:], lev3D[::1,:,:], lev)
        return(v_i[0,:,:])
    elif mode=='th':
        v_i = interpolate(var3D[::1,:,:], lev3D[::1,:,:], lev)
        return(v_i[0,:,:])
    else:
        print('inter2lev only defined for p and th mode.')

def ifs2hgrid(cstfile):
    '''
    IFS analysis/forecast: read horizontal grid from cstfile
    ---
    cstfile: constant file from IFS
    '''
    infile = Dataset(cstfile, mode='r')
    latmin = infile.variables['latmin'][0]
    lonmin = infile.variables['lonmin'][0]
    dellon = infile.variables['dellon'][0]
    dellat = infile.variables['dellat'][0]
    latmax = infile.variables['latmax'][0]+dellon/2
    lonmax = infile.variables['lonmax'][0]+dellat/2
    lon    = np.arange(lonmin,lonmax,dellon)
    lat    = np.arange(latmin,latmax,dellat)
    return lon,lat

def ifs2akbk(cstfile):
    '''
    IFS analysis/forecast: get ak, bk from constant file
    ---
    cstfile: constant file from IFS
    '''
    ak  = readcdf(cstfile,'aklay')
    bk  = readcdf(cstfile,'bklay')
    return ak,bk
    
def era52akbk(pfile,nlev):
    '''
    ERA5: get ak, bk from hyam and hybm
    ---
    nlev: number of vertical levels of  
    '''
    hyam=readcdf(pfile,'hyam')  # 137 levels
    hybm=readcdf(pfile,'hybm')  #   ''
    ak=hyam[hyam.shape[0]-nlev:] # e.g. only 98 levs are used
    bk=hybm[hybm.shape[0]-nlev:] # reduce to 98 levels
    return ak,bk

def echamakbk(ifile):
    '''
    ECHAM-wiso: get ak, bk from hyam and hybm
    ---
    nlev: number of vertical levels of
    '''
    ak=readcdf(ifile,'hyam')  # all levels used
    bk=readcdf(ifile,'hybm')  #   ''
    return ak,bk
