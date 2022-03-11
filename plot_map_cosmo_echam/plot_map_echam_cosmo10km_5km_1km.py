#!/usr/bin/env python
# coding: utf-8
# ----------------------------------------
# mapplots echam-wiso, cosmo 10, 5, 1km
# ----------------------------------------
# Jan 2022 Leonie Villiger
# ----------------------------------------

# ------------------------
# import modules
# ------------------------
import sys
import numpy as np
# plotting
import matplotlib
#matplotlib.use('agg')
import matplotlib.dates as mdate
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
import matplotlib.colors as colorplt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import matplotlib.patches as patches
# language for labels (e.g., dates)
import locale
# mapplots
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
# pylab
import pylab
# netcdf
#import dypy.netcdf as netcdf
import netCDF4
from netCDF4 import Dataset
# time
import datetime as dt
# list files
import glob
#import scipy.ndimage as ndimage
from functions_read_echam_cosmo import read_echam, read_cosmo

def round_to_multiple(number, multiple):
    return multiple * round(number/multiple)

# -------------------------------
# setting

# !!! parameters !!!!!!!!!!!!! 

# get from bash script ($date $var $lev $domain)
arg      = sys.argv 
date     = arg[1]
timestep = dt.datetime(int(date[0:4]),int(date[4:6]),int(date[6:8]),int(date[-2:]))
var      = arg[2]
lev      = float(arg[3]) #hPa
domain   = arg[4]

## set manually (for testing)
#timestep = dt.datetime(2020,2,7,18) #ECHAM available only every 6h
#print(timestep.strftime('%Y%m%d%H'))
#lev=850. #hPa
#var = 'virtpotT'
#domain = 'domain3'

print('Parameters are set to: ', timestep.strftime('%Y%m%d%H'), var, lev, domain)
# !!! parameters !!!!!!!!!!!!!

# map extract:
if domain=='domain1':
    latmin,latmax = -3,29
    lonmin,lonmax = -78,-38
    xloc = np.arange(-70,-25,15)
    yloc = np.arange(0,30,10)
elif domain=='domain2':
    latmin,latmax = 11,16
    lonmin,lonmax = -61,-55
    xloc = np.arange(-61,-53,2)
    yloc = np.arange(11,17,2)
elif domain=='domain3':
    latmin,latmax = 10.5,16
    lonmin,lonmax = -61,-54.5
    xloc = np.arange(-61,-53,2)
    yloc = np.arange(11,17,2)
    #    latmin,latmax = 12,14.5
#    lonmin,lonmax = -60.5,-57
#    xloc = np.arange(-60,-56,1)
#    yloc = np.arange(12,14.5,0.5)

# path
idir_echam = '../echam_data/' #/cfsdYYYYMMDDHH.nc #2019 and 2020
idir_cos10 = '../cos10_data/' #lffdYYYYMMDDHH.ncfile* [1,2,3,4,5]
idir_cos05 = '../cos05_data/'
idir_cos01 = '../cos01_data/'
atr_file   = '../data_ATR/ATR_allflights_lon_lat_alt_p.nc'
odir       = '../plots_cosmo_nesting/'

# -------------------------------
# read data
print('\nreading data')
data_echam   = read_echam(idir_echam, var, timestep, lev)
data_echam_qc= read_echam(idir_echam, 'QC', timestep, lev) #scaled by 1000.

data_cos10 = read_cosmo(idir_cos10, var, timestep, lev)
data_cos10_qc = read_cosmo(idir_cos10, 'QC', timestep, lev) #scaled by 1000.
data_cos10_qr = read_cosmo(idir_cos10, 'QR', timestep, lev) #scaled by 1000.

data_cos05 = read_cosmo(idir_cos05, var, timestep, lev)
data_cos05_qc = read_cosmo(idir_cos05, 'QC', timestep, lev) #scaled by 1000.
data_cos05_qr = read_cosmo(idir_cos05, 'QR', timestep, lev) #scaled by 1000.

data_cos01 = read_cosmo(idir_cos01, var, timestep, lev)
data_cos01_qc = read_cosmo(idir_cos01, 'QC', timestep, lev) #scaled by 1000.
data_cos01_qr = read_cosmo(idir_cos01, 'QR', timestep, lev) #scaled by 1000.

# find min and max of all data to set cmap level
if domain=='domain1':
    qu = 0.99 #upper quantile
    ql = 0.01 #lower q.
elif domain in ['domain2','domain3']:
    qu = 1.
    ql = 0.

idx_echam = np.where((data_echam['lon']>=lonmin)&(data_echam['lon']<=lonmax)&(data_echam['lat']>=latmin)&(data_echam['lat']<=latmax))
min_echam = np.nanquantile(data_echam[var][idx_echam[0],idx_echam[1]], q=ql)
max_echam = np.nanquantile(data_echam[var][idx_echam[0],idx_echam[1]], q=qu)
del idx_echam

idx_cos10 = np.where((data_cos10['lon']>=lonmin)&(data_cos10['lon']<=lonmax)&(data_cos10['lat']>=latmin)&(data_cos10['lat']<=latmax))
min_cos10 = np.nanquantile(data_cos10[var][idx_cos10[0],idx_cos10[1]], q=ql)
max_cos10 = np.nanquantile(data_cos10[var][idx_cos10[0],idx_cos10[1]], q=qu)
del idx_cos10

idx_cos05 = np.where((data_cos05['lon']>=lonmin)&(data_cos05['lon']<=lonmax)&(data_cos05['lat']>=latmin)&(data_cos05['lat']<=latmax))
min_cos05 = np.nanquantile(data_cos05[var][idx_cos05[0],idx_cos05[1]], q=ql)
max_cos05 = np.nanquantile(data_cos05[var][idx_cos05[0],idx_cos05[1]], q=qu)
del idx_cos05

idx_cos01 = np.where((data_cos01['lon']>=lonmin)&(data_cos01['lon']<=lonmax)&(data_cos01['lat']>=latmin)&(data_cos01['lat']<=latmax))
min_cos01 = np.nanquantile(data_cos01[var][idx_cos01[0],idx_cos01[1]], q=ql)
max_cos01 = np.nanquantile(data_cos01[var][idx_cos01[0],idx_cos01[1]], q=qu)
del idx_cos01

# prepare levels for colormap
vmin   = np.ceil(np.nanmin((min_echam, min_cos10, min_cos05, min_cos01)))
vmax   = np.floor(np.nanmax((max_echam, max_cos10, max_cos05, max_cos01)))
#vmin   = round_to_multiple(vmin, 2)
#vmax   = round_to_multiple(vmax, 2)
levels = np.linspace(vmin,vmax,50,endpoint=True)
print('levels: ', levels[0], levels[-1])

# read lon,lat of all ATR flights
ncfile  = Dataset(atr_file)
atr_lon = np.array(ncfile.variables['lon'])
atr_lat = np.array(ncfile.variables['lat'])
del ncfile

# -------------------------------
# general plotting settings
fs=12 #fontsize
latbco, lonbco = 13.16, -59.43

pltdict={'QV'      : ['YlGnBu', 'both', 'Specific humidity', '[g kg$^{-1}$]'],
         'dD'      : ['Spectral', 'both', '$\delta^{2}$H$_\mathrm{vapour}$', '[‰]'],
         'd18O'    : ['Spectral', 'both', '$\delta^{18}$O$_\mathrm{vapour}$', '[‰]'],
         'dexc'    : ['Spectral_r', 'both', 'd-excess$_\mathrm{vapour}$', '[‰]'],
         'RH'      : ['Blues', 'both', 'Relative humidity', '[%]'],
         'T'       : ['RdYlBu_r', 'both', 'Temperature', '[°C]'],
         'virtpotT': ['BrBG', 'both', 'Virtual potential temperature', '[°C]']}

# -------------------------------
# plotting
print('\nplotting')

nrow,ncol = 2,2

# prepare map features land and islands
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                    edgecolor='black',
                                    facecolor='none')

# plot data
plt.close('all')
if domain == 'domain1':
    fig = plt.figure(figsize=(ncol*5,nrow*4.5))
    gs  = gridspec.GridSpec(nrows=nrow+1, ncols=ncol, height_ratios=[10, 10, 1])
    plt.subplots_adjust(top=0.964,bottom=0.103,left=0.067,right=0.974,hspace=0.175,wspace=0.054)
elif domain in 'domain2':
    fig = plt.figure(figsize=(ncol*5,nrow*4.8))
    gs  = gridspec.GridSpec(nrows=nrow+1, ncols=ncol, height_ratios=[10, 10, 1])
    plt.subplots_adjust(top=0.97,bottom=0.095,left=0.067,right=0.962,hspace=0.144,wspace=0.147)
elif domain in 'domain3':
    fig = plt.figure(figsize=(ncol*5,nrow*4.2))
    gs  = gridspec.GridSpec(nrows=nrow+1, ncols=ncol, height_ratios=[10, 10, 1])

# distribute axes
ax1  = fig.add_subplot(gs[0,0],projection=ccrs.PlateCarree())
ax2  = fig.add_subplot(gs[0,1],projection=ccrs.PlateCarree())
ax3  = fig.add_subplot(gs[1,0],projection=ccrs.PlateCarree())
ax4  = fig.add_subplot(gs[1,1],projection=ccrs.PlateCarree())
cax  = fig.add_subplot(gs[2,:]) #colorbar
cax.axis('off')

# plot echam
ax1.set_title('(a)', loc='left', fontsize=fs)
ax1.set_title('ECHAMwiso', loc='right', fontsize=fs)
p1 = ax1.contourf(data_echam['lon'], data_echam['lat'], data_echam[var], levels=levels, cmap=pltdict[var][0], transform=ccrs.PlateCarree(),extend=pltdict[var][1])
#p1 = ax1.contourf(data_echam['lon'], data_echam['lat'], data_echam[var], levels = levdict[var][str(int(lev))], cmap=pltdict[var][0], transform=ccrs.PlateCarree(),extend=pltdict[var][1])
ax1.contour(data_echam_qc['lon'], data_echam_qc['lat'], data_echam_qc['QC'], colors='blue', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())

# plot cosmo 10km
ax2.set_title('(b)', loc='left', fontsize=fs)
ax2.set_title('COSMO$_\mathrm{iso,10km}$', loc='right', fontsize=fs)
ax2.contourf(data_cos10['lon'], data_cos10['lat'], data_cos10[var], levels=levels,
    cmap=pltdict[var][0], transform=ccrs.PlateCarree(),extend=pltdict[var][1])
ax2.contour(data_cos10_qc['lon'], data_cos10_qc['lat'], data_cos10_qc['QC'], colors='blue', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())
ax2.contour(data_cos10_qr['lon'], data_cos10_qr['lat'], data_cos10_qr['QR'], colors='fuchsia', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())

# plot cosmo 05km
ax3.set_title('(c)', loc='left', fontsize=fs)
ax3.set_title('COSMO$_\mathrm{iso,5km}$', loc='right', fontsize=fs)
ax3.contourf(data_cos05['lon'], data_cos05['lat'], data_cos05[var], levels=levels,
    cmap=pltdict[var][0], transform=ccrs.PlateCarree(),extend=pltdict[var][1])
ax3.contour(data_cos05_qc['lon'], data_cos05_qc['lat'], data_cos05_qc['QC'], colors='blue', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())
ax3.contour(data_cos05_qr['lon'], data_cos05_qr['lat'], data_cos05_qr['QR'], colors='fuchsia', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())

# plot cosmo 01km
ax4.set_title('(d)', loc='left', fontsize=fs)
ax4.set_title('COSMO$_\mathrm{iso,1km}$', loc='right', fontsize=fs)
ax4.contourf(data_cos01['lon'], data_cos01['lat'], data_cos01[var], levels=levels,
    cmap=pltdict[var][0], transform=ccrs.PlateCarree(),extend=pltdict[var][1])
ax4.contour(data_cos01_qc['lon'], data_cos01_qc['lat'], data_cos01_qc['QC'], colors='blue', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())
ax4.contour(data_cos01_qr['lon'], data_cos01_qr['lat'], data_cos01_qr['QR'], colors='fuchsia', levels=[0.001],
        linewidths=1, linestyles='-', transform=ccrs.PlateCarree())

# settings true for all plot axes
for ax in [ax1, ax2, ax3, ax4]:
    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree()) # for selection region
    
    # grid ticks and lines
    ax.set_xticks(xloc, crs=ccrs.PlateCarree());
    ax.set_yticks(yloc, crs=ccrs.PlateCarree());
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelsize=fs)
#    ax.gridlines(color='gray',alpha=0.5, linestyle='--',linewidth=1.5)
    if ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels(), visible=False)
    if ax in [ax2,ax4]:
        plt.setp(ax.get_yticklabels(), visible=False)

    if domain=='domain1':
        # add land/country outline
        ax.add_feature(land_50m)
        
#        # location BCO
#        ax.plot(lonbco, latbco, marker='o', markerfacecolor='none', markeredgecolor='k', markersize=15,
#            linewidth=1, transform=ccrs.PlateCarree(),label='BCO')
    
    elif domain in ['domain2','domain3']:
        # add land/country outline
        ax.add_feature(land_50m, facecolor='lightgray', edgecolor='darkgray')
    
#        # atr flight path
#        ax.scatter(atr_lon, atr_lat, marker='.', color='k', alpha=0.1)
    
    # add square embracing atr flight path
    square1 = patches.Rectangle((np.nanmin(atr_lon), np.nanmin(atr_lat)),
        np.nanmax(atr_lon)-np.nanmin(atr_lon), np.nanmax(atr_lat)-np.nanmin(atr_lat),
        linestyle='-', linewidth=1, edgecolor='red', facecolor='none')
    ax.add_patch(square1)


# colorbar
ticks = [vmin,vmin+(vmax-vmin)*0.25, vmin+(vmax-vmin)*0.5, vmin+(vmax-vmin)*0.75, vmax]
cbar = plt.colorbar(p1, ax=cax, orientation='horizontal',fraction=0.9, pad=0.07, shrink=0.7, ticks=ticks)
cbar.ax.tick_params(labelsize=fs)
#plt.setp(cbar.ax.get_xticklabels()[::2], visible=False) #set every second tick invisible
cbar.set_label(pltdict[var][2]+' at '+str(int(lev))+'hPa '+pltdict[var][3]+' on '+timestep.strftime('%d %b %Y, %H UTC'),fontsize=fs)

# save the figure
#plt.tight_layout()
plt.savefig(odir+'map_echam_cosmo_'+var+'_'+str(int(lev))+'hPa_'+timestep.strftime('%Y%m%d_%H')+'_'+domain+'.png', bbox_inches='tight')
plt.show()
