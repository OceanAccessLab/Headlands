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
import hl_modules as hl  
from matplotlib.colors import from_levels_and_colors
import unicodedata
import os

def is_number(s):
    '''
    Used for differentiate numbers from letters in scorecards.
    https://www.pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
    
    '''
    try:
        float(s)
        return True
    except ValueError:
        pass 
    try:
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass 
    return False

## ---- Headlands scorecard ---- ##
df = hl.annual_std_anom('Arnolds_Cove.nc')
df = pd.concat([df, hl.annual_std_anom('Bristols_Hope.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Comfort_Cove.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Hampden.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Lumsden_5m.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Lumsden_15m.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Melrose_5m.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Melrose_15m.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Old_Bonaventure.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Stock_Cove.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Upper_Gullies.nc')], axis=1, sort=False)
df = pd.concat([df, hl.annual_std_anom('Winterton.nc')], axis=1, sort=False)

df.columns = ["Arnold's Cove",
    "Bristol's Hope",
    "Comfort Cove",
    "Hampden",
    "Lumsden 5m",
    "Lumsden 15m",
    "Melrose 5m",
    "Melrose 15m",
    "Old Bonaventure",
    "Stock Cove",
    "Upper Gullies",
    "Winterton"]

# Save .csv for archive
df.to_csv('Headlands_std_anom.csv', sep=',', float_format='%0.3f')

# ignore MEAN/SD for now
std_anom = df.iloc[0:-2]
year_list = std_anom.index.astype('str')
year_list = [i[2:4] for i in year_list] # 2-digit year

# Transpose for table
std_anom = std_anom.T
std_anom['MEAN'] = df.loc['MEAN']
std_anom['SD'] = df.loc['SD']
std_anom.rename(columns={'MEAN': r'$\rm \overline{x}$', 'SD': r'sd'}, inplace=True)

# Add Headlands index
std_anom.at['Headlands index'] = std_anom.mean() 
std_anom.iloc[-1,-1]=np.nan
std_anom.iloc[-1,-2]=np.nan  

# Get text values +  cell color
year_list.append(r'$\rm \overline{x}$') # add 2 extra columns
year_list.append(r'sd')   
vals = np.around(std_anom.values,1)
vals[vals==-0.] = 0.
vals_color = vals.copy()
vals_color[:,-1] = 0 # No color to last two columns (mean and STD)
vals_color[:,-2] = 0

#vals_color[(vals_color<0.5) & (vals_color>-.5)] = 0.

# Build the colormap
vmin = -3.49
vmax = 3.49
midpoint = 0
levels = np.linspace(vmin, vmax, 15)
midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
normal = plt.Normalize(-3.49, 3.49)
reds = plt.cm.Reds(np.linspace(0,1, num=7))
blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
whites = [(1,1,1,1)]*2
colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
cmap, norm = from_levels_and_colors(levels, colors, extend='both')
cmap_r, norm_r = from_levels_and_colors(levels, np.flipud(colors), extend='both')

nrows, ncols = std_anom.index.size+1, std_anom.columns.size
hcell, wcell = 0.5, 0.5
hpad, wpad = 1, 1    

# Here figure
fig=plt.figure(figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
ax = fig.add_subplot(111)
ax.axis('off')
#do the table
header = ax.table(cellText=[['']],
                      colLabels=['--  Headlands Coastal Temperature Network (May-Oct)  --'],
                      loc='center'
                      )
header.set_fontsize(13)
the_table=ax.table(cellText=vals, rowLabels=std_anom.index, colLabels=year_list,
                    loc='center', cellColours=cmap(normal(vals_color)), cellLoc='center',
                    bbox=[0, 0, 1, 0.5]
                    )
# change font color to white where needed:
the_table.auto_set_font_size(False)
the_table.set_fontsize(12.5)
table_props = the_table.properties()
table_cells = table_props['child_artists']
last_columns = np.arange(vals.shape[1]-2, vals.shape[1]) # last columns
for key, cell in the_table.get_celld().items():
    cell_text = cell.get_text().get_text()
    if is_number(cell_text) == False:
        pass
    elif key[0] == 0: #year's row = no color
        pass
    elif (cell_text=='nan'):
        cell._set_facecolor('lightgray')
        cell._text.set_color('lightgray')
    elif key[1] in last_columns:
         cell._text.set_color('darkslategray')
    elif (np.float(cell_text) <= -1.5) | (np.float(cell_text) >= 1.5) :
        cell._text.set_color('white')

    # index in bold
    if key[0] == 13:
        cell._text.set_fontweight('bold')

filename = 'scorecards_Headlands.png'
plt.savefig(filename, dpi=300)
os.system('convert -trim ' + filename + ' ' + filename)
