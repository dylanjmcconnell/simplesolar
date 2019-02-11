import numpy as np
import pandas as pd
import xarray as xr
import os
import datetime
import time
import progressbar

def get_npdatetime(pathname):
    """Accepts a pathname (including file and extension) as a string and returns the datetime value associated with UTS timezone. Assumed to be in format: 'solar_dni_20150101_01_UT.txt'"""
    filename = pathname[-27:]
    year = filename[10:14]
    month = filename[14:16]
    day = filename[16:18]
    hour = filename[19:21]
    
    return (np.datetime64('{}-{}-{}T{}:00:00'.format(year, month, day, hour)))

def get_datetime(pathname):
    """Accepts a filename as a string and returns the datetime value associated with UTS timezone. Assumed to be in format: 'solar_dni_20150101_01_UT.txt'"""
    filename = pathname[-27:]
    year = filename[10:14]
    month = filename[14:16]
    day = filename[16:18]
    hour = filename[19:21]
    
    return (pd.to_datetime('{}-{}-{}T{}:00:00'.format(year, month, day, hour)))


def get_files(mypath, firstnine):
    """Returns a list of files that have filenames starting with "firstnine" as list of lists [fullpath, filename]."""
    fileinfo = []
    print("Getting file list.")
    for root, dirs, files in progressbar.progressbar(os.walk(mypath, topdown = False)):
        for name in files:
            if name[0:9] == firstnine:
                fileinfo.append((get_datetime(name),os.path.join(root,name)))
                
    print (str(datetime.datetime.now()), ": Number of files =", len(fileinfo))
    
    fileinfo.sort()
    return (fileinfo[:])


def create_xarray_list(filelist):
    """Creates xarrays from files in the filelist within the current directory."""
    
    #Create lat/long np.arrays for creation of xarrays
    lats = np.flip(np.arange(-43.95, -10.05, 0.05), axis = 0)
    lons = np.arange(112.05, 153.96, 0.05)
    
    
    print(str(datetime.datetime.now()), ": Creating xarray of files.")
    xarray_list = []
    
    for file in progressbar.progressbar(filelist):
        
        #creates a timestamp
        filedatetime = file[0]
        
        #Creates a np array with additional empty dimension (for time) from the file with appropriate filename.
        np_rad = np.expand_dims(np.loadtxt(file[1], delimiter = ' ', skiprows = 6), axis = 0)
        
        
        #Creates a xarray DataArray with 3 dims, coords time, latitude, longitude.
        xarray_list.append(xr.DataArray(np_rad, coords = [[filedatetime], lats, lons], dims=['time','latitude', 'longitude']))
    
    
    radArray = xr.concat((xarray_list), 'time')
    radArray.attrs['unit'] = 'W/m2'
        
    return(radArray)


def get_months(filelist):
    files = []
    month = []
    n = 0
    for x in filelist:
        n+=1
        month.append(x)
        if n == len(filelist):
            files.append(month)
            return(files)
        
        elif x[0].month != filelist[n][0].month:
            files.append(month)
            month = []


def convert_bomdata(mypath = '/Users/felixsilberstein/Desktop/SOLAR_DATA/', targetpath = '/Users/felixsilberstein/Desktop/sandbox/', targetname = 'radiation'):
    """Saves a netCDF to the path specified."""
    DNI_files = get_months(get_files(mypath, 'solar_dni'))
    GHI_files = get_months(get_files(mypath, 'solar_ghi'))
    
    for x,y in zip(DNI_files, GHI_files):
        date = x[0][0].strftime('%Y%m')
        dniArray = create_xarray_list(x)
        ghiArray = create_xarray_list(y)
    
        radDataset = xr.Dataset({'dni': dniArray, 'ghi' : ghiArray})
        radDataset.attrs['name'] = 'solarradiation'
    
        print(str(datetime.datetime.now()), ": Saving {} as netCDF.".format(date))
        radDataset.to_netcdf('{0}{1}{2}.nc'.format(targetpath, targetname, date), unlimited_dims= 'time')
    
    return(None)






