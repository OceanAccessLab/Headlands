'''
This script will create a netCDF file containing raw deployment from the Headlands project.
Historical data are from trimmed NAFC's native format *.rpf.
Newer data since deployment year 2013 are found in .csv format.

At the time of writing this script,

Ryan.Doody@dfo-mpo.gc.ca was responsible for the Headlands program
and
Charlie.Bishop@dfo-mpo.gc.ca generated historical archive.
Also:
Giulia Bronzi has been working on this during a summer intership (June-August 2020)
and
Jonathan.Coyne@dfo-mpo.gc.ca will take over the annual data processing in 2025.


Frederic.Cyr@dfo-mpo.gc.ca
January/February 2022


'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import xarray as xr
import hl_modules as hl  
import os 

## Some parameters
YEAR = 2024

## Annual daily plot with climatology:
hl.plotClimatology('Arnolds_Cove.nc', YEAR, siteName="Arnold's Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Bristols_Hope.nc', YEAR, siteName="Bristol's Hope", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Comfort_Cove.nc', YEAR, siteName="Comfort Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Hampden.nc', YEAR, siteName="Hampden", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Lumsden_15m.nc', YEAR, siteName="Lumsden (15m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Lumsden_5m.nc', YEAR, siteName="Lumsden (5m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Melrose_15m.nc', YEAR, siteName="Melrose (15m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Melrose_5m.nc', YEAR, siteName="Melrose (5m)", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Old_Bonaventure.nc', YEAR, siteName="Old Bonaventure", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Stock_Cove.nc', YEAR, siteName="Stock Cove", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Upper_Gullies.nc', YEAR, siteName="Upper Gullies", method = 'doy', XLIM=[100,350])
hl.plotClimatology('Winterton.nc', YEAR, siteName="Winterton", method = 'doy', XLIM=[100,350])
plt.close()

# Convert to a subplot
os.system('montage \
ArnoldsCove_' + str(YEAR) + '.png  \
BristolsHope_' + str(YEAR) + '.png \
ComfortCove_' + str(YEAR) + '.png \
Hampden_' + str(YEAR) + '.png \
Lumsden_5m_' + str(YEAR) + '.png \
Lumsden_15m_' + str(YEAR) + '.png \
Melrose_5m_' + str(YEAR) + '.png \
Melrose_15m_' + str(YEAR) + '.png \
OldBonaventure_' + str(YEAR) + '.png \
StockCove_' + str(YEAR) + '.png \
UpperGullies_' + str(YEAR) + '.png \
Winterton_' + str(YEAR) + '.png \
-tile 3x4 -geometry +10+10  -background white  headlands_' + str(YEAR) + '_montage.png')


# Salmon mass die off(doy-240+)
#hl.plotClimatology('Arnolds_Cove.nc', 2019, siteName="Arnold's Cove", method = 'doy')


## Annual monthly barplot with climatology and anomaly:
hl.annual_plot_anom('Arnolds_Cove.nc', "Arnold's Cove", YEAR)
hl.annual_plot_anom('Bristols_Hope.nc', "Bristol's Hope", YEAR)
hl.annual_plot_anom('Comfort_Cove.nc', "Comfort Cove", YEAR)
hl.annual_plot_anom('Hampden.nc', "Hampden", YEAR)
hl.annual_plot_anom('Lumsden_15m.nc', "Lumsden (15m)", YEAR)
hl.annual_plot_anom('Lumsden_5m.nc', "Lumsden (5m)", YEAR)
hl.annual_plot_anom('Melrose_15m.nc', "Melrose (15m)", YEAR)
hl.annual_plot_anom('Melrose_5m.nc', "Melrose (5m)", YEAR)
hl.annual_plot_anom('Old_Bonaventure.nc', "Old Bonaventure", YEAR)
hl.annual_plot_anom('Stock_Cove.nc', "Stock Cove", YEAR)
hl.annual_plot_anom('Upper_Gullies.nc', "Upper Gullies", YEAR)
hl.annual_plot_anom('Winterton.nc', "Winterton", YEAR)
plt.close()

# Convert to a subplot
os.system('montage Arnolds_Cove_monthly.png  Bristols_Hope_monthly.png Comfort_Cove_monthly.png Hampden_monthly.png Lumsden_5m_monthly.png Lumsden_15m_monthly.png Melrose_5m_monthly.png Melrose_15m_monthly.png Old_Bonaventure_monthly.png Stock_Cove_monthly.png Upper_Gullies_monthly.png Winterton_monthly.png -tile 3x4 -geometry +10+10  -background white  headlands_monthly_montage.png')


## Scorecards figure
%my_run hl_scorecards.py
