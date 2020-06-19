'''
AZMP script to read Headlands thermograph data from D. Senciall's .rpf format

data in : /home/cyrf0006/data/Headlands/
process in : /home/cyrf0006/AZMP/Headlands

Frederic.Cyr@dfo-mpo.gc.ca - February 2019

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import getpass
import xarray
import hl_modules

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

    
#makes list of every file in ComfortCove folder
fileList = hl_modules.getFileList("ComfortCove")

#reads the rpf files listed in fileList
df_all = hl_modules.readrpf(fileList)

#removes duplicates from df_all
df_all = hl_modules.removeDuplicates(df_all)

headerInfo = hl_modules.extractHeaders(fileList)



# monthly average
df_monthly = df_all.resample('M').mean()
df_monthly.plot()
plt.title("Monthly Average")
plt.show()

df_monthly.to_csv('comfort_cove_thermograph_1989-2017_monthly.csv')


# #June-July only:
# df_summer = df_monthly[(df_monthly.index.month>=6) & (df_monthly.index.month<=7)]
# df_summer = df_summer.resample('As').mean()
# df_summer.index = df_summer.index.year
# df_summer.to_csv('comfort_cove_thermograph_1989-2017_June-July.csv')


# ## ---- plot summer data in anomalies ---- ##
# anom = (df_summer - df_summer.mean()) / df_summer.std()
# df1 = anom[anom<0]
# df2 = anom[anom>0]
# fig = plt.figure(4)
# fig.clf()
# width = .9
# p1 = plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, color='steelblue')
# p2 = plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, alpha=0.8, color='indianred')
# plt.ylabel('Standardized Anomaly')
# #plt.xlabel('Year')
# plt.title('Comfort Cove temperature (June-July)')
# plt.grid()
# fig.set_size_inches(w=15,h=9)
# fig_name = 'Comfort_Cove_anomalies.png'
# #plt.annotate('data source: NCDC/NOAA', xy=(.75, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
# fig.savefig(fig_name, dpi=300)
# os.system('convert -trim ' + fig_name + ' ' + fig_name)


## ---- annual curve ---- ##
#originally looking at 2018 but wasnt in file, using 2011 instead
df_all['woy'] = df_all.index.weekofyear
weekly_clim = df_all.groupby('woy').mean()
weekly_std = df_all.groupby('woy').std()
df_2011 = df_all[df_all.index.year>=2011]
weekly_2011 = df_2011.groupby('woy').mean()

ax = weekly_clim.plot(linewidth=3, legend=None)
weekly_2011.plot(ax=ax, linewidth=3)
plt.fill_between(weekly_clim.index, 
                  np.squeeze(weekly_clim.values+weekly_std.values),
                  np.squeeze(weekly_clim.values-weekly_std.values), 
                  facecolor='steelblue', 
                  interpolate=True , 
                  alpha=.3)
plt.ylabel(r'T ($^{\circ}C$)')
plt.xlabel('Week of the year')
plt.title('Comfort Cove temperature')
plt.xlim([0,53])
plt.ylim([-2,18])
plt.grid()
plt.legend(['1989-2018 average', '2011']) #not correct years
fig = ax.get_figure()
fig.set_size_inches(w=12,h=8)
fig_name = 'Comfort_Cove_T.png'
fig.savefig(fig_name, dpi=200)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



