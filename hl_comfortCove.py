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
import pandas as pd
import os

# Adjust fontsize/weight
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

# Generate the list
infiles = 'comfortcove.list'
os.system('ls ~/Desktop/ComfortCove/*.rpf > ' + infiles)
filelist = np.genfromtxt(infiles, dtype=str)
filelist = np.reshape(filelist, filelist.size) 

dfs = []
for fname in filelist: 
    df = pd.read_csv(fname, sep='\s+',  parse_dates={'datetime': [0, 1]}, header=15)
    df = df.set_index('datetime')
    df.columns = ['temperature']
    df = df.replace(9999.99, np.NaN)
    #print(fname, df.max())
    dfs.append(df)
    
# concatenate all data <----------- Here we need to check if duplicates   
df_all = pd.concat(dfs, axis=0)
df_all = df_all.sort_index()

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
plt.fill_between(weekly_clim.index, np.squeeze(weekly_clim.values+weekly_std.values),np.squeeze(weekly_clim.values-weekly_std.values), facecolor='steelblue', interpolate=True , alpha=.3)
plt.ylabel(r'T ($^{\circ}C$)')
plt.xlabel('Week of the year')
plt.title('Comfort Cove temperature')
plt.xlim([0,53])
plt.ylim([-2,18])
plt.grid()
plt.legend(['1989-2018 average', '2011']) #not correct years
fig.set_size_inches(w=15,h=9)
fig_name = 'Comfort_Cove_T.png'
fig.savefig(fig_name, dpi=300)
os.system('convert -trim ' + fig_name + ' ' + fig_name)



dfsHeader = [] 

#Creates an array of dataframes, one for each file. Each dataframe is composed of the file's header.
for fname in filelist:
    df = pd.read_csv(fname, nrows = 14)
    headerValue = df["HEADER"].str.split("=", expand = True)
    df["Title"] = headerValue[0]
    df["Value"] = headerValue[1]
    df = df[['Title',"Value"]]
    dfsHeader.append(df)

print(dfsHeader[1]) #first file's header

 
#array for each header component   
station = []
siteName = []
startDate = []
startTime = []
endDate = []
endTime = []
latitude = []
longitude = []
instType = []
serialNumber = []
waterDepth = []
instDepth = []
samplingInterval = []
fileName = []


#extracts each individual value from each header item, from each file
#adds each value to correlated array
count = 0
for file in dfsHeader:
    station.append(dfsHeader[count]["Value"][0])
    siteName.append(dfsHeader[count]["Value"][1])
    startDate.append(dfsHeader[count]["Value"][2])
    startTime.append(dfsHeader[count]["Value"][3])
    endDate.append(dfsHeader[count]["Value"][4])
    endTime.append(dfsHeader[count]["Value"][5])
    latitude.append(dfsHeader[count]["Value"][6])
    longitude.append(dfsHeader[count]["Value"][7])
    instType.append(dfsHeader[count]["Value"][8])
    serialNumber.append(dfsHeader[count]["Value"][9])
    waterDepth.append(dfsHeader[count]["Value"][10])
    instDepth.append(dfsHeader[count]["Value"][11])
    samplingInterval.append(dfsHeader[count]["Value"][12])
    fileName.append(dfsHeader[count]["Value"][13])
    count += 1
    
    

headers = {'Station': station,
           'Site Name': siteName,
           'Start Date': startDate,
           'Start Time': startTime,
           'End Date': endDate,
           'End Time': endTime,
           'Latitude': latitude,
           'Longitude': longitude,
           'Inst Type': instType,
           'Serial Number': serialNumber,
           'Water Depth': waterDepth,
           'Inst Depth': instDepth,
           'Sampling Interval': samplingInterval,
           'File Name': fileName }

#creates new dataframe based on the arrays
headersdf = pd.DataFrame(headers)

#dataframe where each column is a component of the header
#each row is an individual file
print(headersdf[:5]) #first five rows