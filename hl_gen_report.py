'''
This script will create a netCDF file containing raw deployment from the Headlands project.
Historical data are from trimmed NAFC's native format *.rpf.
Newer data since deployment year 2013 are found in .csv format.

At the time of writing this script,

Ryan.Doody@dfo-mpo.gc.ca was responsible for the Headlands program
and
Charlie.Bishop@dfo-mpo.gc.ca generated historical archive.

This script was largely written by Giulia Bronzi in June-August 2020


Frederic.Cyr@dfo-mpo.gc.ca
January/February 2022


'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
## import datetime
## import os
## import xarray
## import copy
## import seaborn as sns
## from datetime import date
## from scipy.stats import linregress
## from dateutil.relativedelta import relativedelta

import xarray as xr
import hl_modules as hl  
import os 
## ## Input params to be put in function
## site = 'Arnolds_Cove'; path_name=None
## site = 'Bristols_Hope'; path_name=None
## site = 'ComfortCove'; path_name=None
## site = 'Hampden'; path_name=None
## site = 'Lumsden'; path_name='/5m'
## site = 'Lumsden'; path_name='/15m'
## site = 'Melrose'; path_name='/5m'
## site = 'Melrose'; path_name='/15m'
## site = 'OldBonaventure'; path_name=None
## site = 'StockCove'; path_name=None
## site = 'UpperGullies'; path_name=None
## site = 'Winterton'; path_name=None
## ## 
## filelist = hl.getFileList(site, path_name=path_name)

## # Read .rpf files
## filelist = hl.getFileList('Arnolds_Cove')
## df_all = hl.readrpf(filelist)

## Quick plot:
#ds = xr.open_dataset('Arnolds_Cove.nc')   
#df = ds.temperature.to_pandas()
#df.plot()

## Annual daily plot with climatology:
hl.plotClimatology('Arnolds_Cove.nc', 2021, siteName="Arnold's Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Bristols_Hope.nc', 2021, siteName="Bristol's Hope", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Comfort_Cove.nc', 2021, siteName="Comfort Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Hampden.nc', 2021, siteName="Hampden", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Lumsden_15m.nc', 2021, siteName="Lumsden (15m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Lumsden_5m.nc', 2021, siteName="Lumsden (5m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Melrose_15m.nc', 2021, siteName="Melrose (15m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Melrose_5m.nc', 2021, siteName="Melrose (5m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Old_Bonaventure.nc', 2021, siteName="Old Bonaventure", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Stock_Cove.nc', 2021, siteName="Stock Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Upper_Gullies.nc', 2021, siteName="Upper Gullies", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Winterton.nc', 2021, siteName="Winterton", method = 'doy', XLIM=[100,350])

# Convert to a subplot
os.system('montage ArnoldsCove_2021.png  BristolsHope_2021.png ComfortCove_2021.png Hampden_2021.png Lumsden_5m_2021.png Lumsden_15m_2021.png Melrose_5m_2021.png Melrose_15m_2021.png OldBonaventure_2021.png StockCove_2021.png UpperGullies_2021.png Winterton_2021.png -tile 3x4 -geometry +10+10  -background white  headlands_2021_montage.png')


# Salmon mass die off(doy-240+)
hl.plotClimatology('Arnolds_Cove.nc', 2019, siteName="Arnold's Cove", method = 'doy')


## Annual monthly barplot with climatology and anomaly:
hl.annual_plot_anom('Arnolds_Cove.nc', "Arnold's Cove", 2021)
hl.annual_plot_anom('Bristols_Hope.nc', "Bristol's Hope", 2021)
hl.annual_plot_anom('Comfort_Cove.nc', "Comfort Cove", 2021)
hl.annual_plot_anom('Hampden.nc', "Hampden", 2021)
hl.annual_plot_anom('Lumsden_15m.nc', "Lumsden (15m)", 2021)
hl.annual_plot_anom('Lumsden_5m.nc', "Lumsden (5m)", 2021)
hl.annual_plot_anom('Melrose_15m.nc', "Melrose (15m)", 2021)
hl.annual_plot_anom('Melrose_5m.nc', "Melrose (5m)", 2021)
hl.annual_plot_anom('Old_Bonaventure.nc', "Old Bonaventure", 2021)
hl.annual_plot_anom('Stock_Cove.nc', "Stock Cove", 2021)
hl.annual_plot_anom('Upper_Gullies.nc', "Upper Gullies", 2021)
hl.annual_plot_anom('Winterton.nc', "Winterton", 2021)

# Convert to a subplot
os.system('montage Arnolds_Cove_monthly.png  Bristols_Hope_monthly.png Comfort_Cove_monthly.png Hampden_monthly.png Lumsden_5m_monthly.png Lumsden_15m_monthly.png Melrose_5m_monthly.png Melrose_15m_monthly.png Old_Bonaventure_monthly.png Stock_Cove_monthly.png Upper_Gullies_monthly.png Winterton_monthly.png -tile 3x4 -geometry +10+10  -background white  headlands_monthly_montage.png')

