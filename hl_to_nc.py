'''
This script will create a netCDF file containing raw deployment from the Headlands project.
Historical data are from trimmed NAFC's native format *.rpf.
Newer data since deployment year 2013 are found in .csv format.

At the time of writing this script,

Ryan.Doody@dfo-mpo.gc.ca was responsible for the Headlands program
and
Charlie.Bishop@dfo-mpo.gc.ca generated historical archive.

This script was largely written by Giulia Bronzi in June-August 2020

Path example:
/home/cyrf0006/data/Thermograph/DataFiles/Headlands/Melrose/15m/ValidSampleSize/trimmed

To be run in /home/cyrf0006/AZMP/Headlands

usage example to open nc file:
import xarray as xr
ds = xr.open_dataset('Arnolds_Cove.nc')   
df = ds.temperature.to_pandas()
df.plot()


Frederic.Cyr@dfo-mpo.gc.ca
January/February 2022


'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import xarray as xr
import copy
import seaborn as sns
from datetime import date
from scipy.stats import linregress
from dateutil.relativedelta import relativedelta
import netCDF4 as nc
import time as tt

import os
import hl_modules as hl 

## ---- Input params to be put in function ---- ##
#site = 'Arnolds_Cove'; site_tmp = 'Arnold'
site = 'Bristols_Hope'; site_tmp = 'Bristol'
#site = 'Comfort_Cove'; site_tmp = 'Comfort'
#site = 'Hampden'; site_tmp = 'Hampden'
#site = 'Lumsden_5m'; site_tmp = 'Lumsden.*\ 5m'
#site = 'Lumsden_15m'; site_tmp = 'Lumsden.*15m'
#site = 'Melrose_5m'; site_tmp = 'Melrose.*\ 5m'
#site = 'Melrose_15m'; site_tmp = 'Melrose.*15m'
#site = 'Old_Bonaventure'; site_tmp = 'Bonaventure'
#site = 'Stock_Cove'; site_tmp = 'Stock'
#site = 'Upper_Gullies'; site_tmp = 'Foxtrap'
#site = 'Winterton'; site_tmp = 'Winterton'

years_to_append = [2019, 2020, 2021] 


## ---- Read .rpf files + clean ---- #
# Read file list
filelist = hl.getFileList(site, pre_path='~/data/Headlands_Trimmed')
df_all = hl.readrpf(filelist)

## ---- Append more recent files ---- ##
df_update = pd.DataFrame([])
for year in years_to_append:
    path_tmp = r'"/home/cyrf0006/data/Thermographs/' + str(year) + ' Thermographs/Headlands ' + str(year) + '/Trimmed Data"'
    files_tmp = os.popen('grep -l ' + site_tmp + ' ' + os.path.join(path_tmp, 'Minilog*.csv')).read().split('\n')
    #if file_tmp:
    for file_tmp in files_tmp:
        if file_tmp:
            df_tmp = pd.read_csv(os.path.join(path_tmp, file_tmp), header=7, names=['date','time','temperature'], parse_dates={'datetime': ['date', 'time']}, index_col='datetime')
            df_update = df_update.append(df_tmp, sort=True)
    
# Just in case, resample
#df_update = df_update.resample('H').mean()

# update historical data
df_all = df_all.append(df_update)

# Remove duplicates by hourly average
df = df_all.resample('H').mean()


## ---- get meta data and attribute values ---- ##
# From info file:
df_info = pd.read_csv(site + '.info', header=None, sep='=', names=['attribute', 'value'])
df_info.set_index('attribute', inplace=True)
station_name = df_info.loc['station_short_name']
station_long_name = df_info.loc['station_long_name']
station_number = df_info.loc['station_number']
station_latitude = df_info.loc['station_lat']
station_longitude = df_info.loc['station_lon']
instrument_depth = df_info.loc['instrument_depth']

# From file headers (not everything is used)
header = hl.extractHeaders(filelist)
station = str(header['Station'].unique())
siteName = str(header['SiteName'].unique())
start = str((header['StartDate'] + ' ' + header['StartTime']).values)
end = str((header['EndDate'] + ' ' + header['EndTime']).values)
latitude = str(header['Latitude'].unique())
longitude = str(header['Longitude'].unique())
instType = str(header['InstType'].unique())
serialNumber = str(header['SerialNumber'].unique())
waterDepth = str(header['WaterDepth'].unique())
instDepth = str(header['InstDepth'].unique())
samplingInterval = str(header['SamplingInterval'].unique())
fileName = str(header['FileName'].unique())


## ------ Building netCDF file (inspired from MEOPAR & NCAR examples) ------ ##
nc_outfile = site + '.nc'

# File name + global attributes
nc_out = nc.Dataset(nc_outfile, 'w')
nc_out.Conventions = 'CF-1.6'
nc_out.title = 'Historical temperatures at ' + df_info.loc['station_long_name'].value
nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
nc_out.source = 'https://github.com/oceanaccess/Headlands'
nc_out.references = 'https://azmp-nl.github.io/Headlands'
nc_out.description = 'Headlands: Coastal temperature monitoring program in Newfoundland'
nc_out.author = 'Frederic.Cyr@dfo-mpo.gc.ca'
nc.history = 'Created ' + tt.ctime(tt.time())
nc_out.comment = 'Just a trial at the moment, no distribution!'

# Create dimensions
time = nc_out.createDimension('time', None)
level = nc_out.createDimension('level', len(instrument_depth))

# Create coordinate variables
times = nc_out.createVariable('time', np.float64, ('time',))
levels = nc_out.createVariable('level', np.int32, ('level',))

# Create time variable (1-D)
temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

# Create level variables (0-D) | Should this be attributes?
latitudes = nc_out.createVariable('latitude', np.float32, ('level'), zlib=True)
longitudes = nc_out.createVariable('longitude', np.float32, ('level'), zlib=True)
stations = nc_out.createVariable('station_name', np.str, ('level'), zlib=True)
ID = nc_out.createVariable('station_ID', np.int32, ('level'), zlib=True)

# Variable Attributes
latitudes.units = 'degree_north'
longitudes.units = 'degree_east'
times.units = 'hours since 1900-01-01 00:00:00'
times.calendar = 'gregorian'
levels.units = 'm'
levels.standard_name = "depth"
levels.valid_min = 0
temp.units = 'Celsius'
temp.long_name = "Water Temperature" # (may be use to label plots)
temp.standard_name = "sea_water_temperature"

# Fill variables
latitudes[:] = np.array(df_info.loc['station_lat'])
longitudes[:] = np.array(df_info.loc['station_lon'])
stations[:] = np.array(df_info.loc['station_short_name'])
ID[:] = np.array(df_info.loc['station_number'])
temp[:] = df.values

# Fill dimensions
times[:] = nc.date2num(df.index.to_pydatetime(), units = times.units, calendar = times.calendar)
levels[:] = int(instrument_depth.value)

# Close file
nc_out.close()
print('Done!')

