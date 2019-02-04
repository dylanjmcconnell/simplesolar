import numpy as np
import pandas as pd
import xarray as xr
import os
import datetime

import getdata as gd

def get_npdatetime(filename):
    """Accepts a filename as a string and returns the datetime value associated with UTS timezone. Assumed to be in format: 'solar_dni_20150101_01_UT.txt'"""
    
    year = filename[10:14]
    month = filename[14:16]
    day = filename[16:18]
    hour = filename[19:21]
    
    return (np.datetime64('{}-{}-{}T{}:00:00'.format(year, month, day, hour)))

def get_files(mypath, firstnine)
    """Returns a list of files in 'mypath' where the first 9 characters match 'firstnine'."""
    path_list = [f for f in os.listdir(mypath) if f[0:9]==firstnine]
    return(path_list)


def create_xarr(filelist):
    #Create lat/long np.arrays for creation of xarrays
    lats = np.flip(np.arange(-43.95, -10.05, 0.05), axis = 0)
    lons = np.arange(112.05, 153.96, 0.05)

    #Create array for all of the radiation files, also create a datetime list, and 

    for file in filelist:
        #creates a timestamp
        datetime = gd.get_npdatetime(file)
    
        #Creates a np array with additional empty dimension (for time) from the file with appropriate filename.
        np_rad = np.expand_dims(np.loadtxt(file, delimiter = ' ', skiprows = 6), axis = 0)
    
        #Creates a xarray DataArray with 3 dims, coords time, latitude, longitude.
        temp_arr = xr.DataArray(np_rad, coords = [[datetime], lats, lons], dims=['time','latitude', 'longitude'])
    
        try:
            #if the DataArray array exists, add the new temp array to it as additional entry on the time axis.
            radArray = xr.concat((radArray, temp_arr), 'time')
    
        except:
            #if the DataArray doesnt exist then initialise it.
            radArray = temp_arr
        
    radArray = radArray.sortby('time')
    radArray.attrs['unit'] = 'W/m2'
        
    return(radArray)

def combine_xarr(name, dniArray, ghiArray)
    radDataset = xr.Dataset({'dni': dniArray, 'ghi' : ghiArray})
    radDataset.attrs['name'] = name
    return (radDataset)

def save_dataset(dataset)
    radDataset.to_netcdf('radiation.nc', unlimited_dims= 'time')


def convert_bomdata(mypath):
    DNI_files = get_files(mypath, 'solar_dni')
    GHI_files = get_files(mypath, 'solar_ghi')
    
    dniArray = create_xarr(DNI_files)
    ghiArray = create_xarr(GHI_files)
    
    radDataset = xr.Dataset({'dni': dniArray, 'ghi' : ghiArray})
    radDataset.attrs['name'] = 'solarradiation'
    
    radDataset.to_netcdf('radiation.nc', unlimited_dims= 'time')
    
    return(None)



