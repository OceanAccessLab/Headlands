"""
Modules for headlands project:
    1. Get file list from a folder
    2. Read .rpf files
    3. Remove duplicates
    4. Extract header info
    5. Plot annual curve
    6. Plot anomalies
    7. Plot monthly time series
    8. Plot daily time series
    9. Plot climatology 
    10. Plot/find upwells
    12. Upwell table
    12. Plot derivative curve
    13. Write to netCDFk


"""

import os
import getpass
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray
import calendar as cal
import seaborn as sns
import copy


# =============================================================================
# Takes in a folder name and creates a list of each file name.
#
# Returns a list of file names.
# =============================================================================
  
def getFileList(folderName):
    # Check user's path
    if getpass.getuser() == 'cyrf0006':
        files_path = '~/data/Headlands_Trimmed/'
    else:
        files_path = '~/Desktop'
    # Generate the list
    infiles = folderName + '.list'
    os.system('ls ' + os.path.join(files_path, folderName + '/*.rpf') + ' > ' + infiles)
    filelist = np.genfromtxt(infiles, dtype=str)
    filelist = np.reshape(filelist, filelist.size) 
    
    return filelist
 
# =============================================================================
# Takes in a list of file names and extracts the data, creates an array of 
# dataframes (one df for each file) and then concatenates the array into one
# dataframe. 
#
# Returns a dataframe containing all the data from a site.
# =============================================================================
    
def readrpf(filelist):     
    dfs = []
    for fname in filelist: 
        try: 
            df = pd.read_csv(fname, sep='\s+',  
                             parse_dates={'datetime': [0, 1]}, 
                             header=16)
            df = df.set_index('datetime')
            df.columns = ['temperature']
            df = df.replace(9999.99, np.NaN)
            dfs.append(df)
        except:
            print (fname + ' is empty [ignore file]')
            continue
        
    # concatenates all data 
    df_all = pd.concat(dfs, axis=0)
    df_all = df_all.sort_index()
    
    return df_all


# =============================================================================
# Takes in a dataframe, finds average temp between any duplicates,
# removes the duplicates, and returns the dataframe. 
#
# Returns a dataframe.
# =============================================================================
def removeDuplicates(df):
    
    #sets temperature to the average (if not a duplicate, it doesn't affect temp)
    df['temperature'] = ((df.temperature.resample('H').max() + 
                              df.temperature.resample('H').min())/2)

    #makes datetime a column so its easy to delete duplicates
    df = df.reset_index()

    #drops datetime duplicates
    df = df.drop_duplicates(subset='datetime')
    df = df.set_index('datetime')

    return df


# =============================================================================
# Takes in a list of file names and creates a dataframe containting all 
# header information. Each row of the dataframe relates to one file. 
#   
# Returns a dataframe.
# =============================================================================
def extractHeaders(filelist):

    dfsHeader = [] 
    #Creates an array of dataframes, one for each file. 
    #Each dataframe is composed of the file's header.
    for fname in filelist:
        dfH = pd.read_csv(fname, nrows = 14)
        headerValue = dfH["HEADER"].str.split("=", expand = True)
        dfH["Title"] = headerValue[0]
        dfH["Value"] = headerValue[1]
        dfH = dfH[['Title',"Value"]]
        dfsHeader.append(dfH)
    
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
    
    colNames = ['Station', 'Site Name', 'Start Date', 'Start Time', 'End Date', 
                'End Time','Latitude', 'Longitude',  'Inst Type', 'Serial Number',
                'Water Depth','Inst Depth', 'Sampling Interval', 'File Name']
    
    #creates new dataframe, each column is a component of the header
    #each row represents an individual file
    headersdf = pd.DataFrame(headers, columns = colNames)

    return headersdf


# =============================================================================
# Takes in a dataframe and a year and produces a plot comparing the inputted
# year with the average from 1989-2018.
#
# Plots graph. 
# =============================================================================
def plotAnnualCurve(df, year, siteName):
    
    dfA = df
    
    dfA['woy'] = dfA.index.weekofyear
    weekly_clim = dfA.groupby('woy').mean()
    weekly_std = dfA.groupby('woy').std()
    df_year = dfA[dfA.index.year >= year]
    weekly_year = df_year.groupby('woy').mean()
    
    ax = weekly_clim.plot(linewidth=3, legend=None)
    weekly_year.plot(ax=ax, linewidth=3)
    plt.fill_between(weekly_clim.index, 
                      np.squeeze(weekly_clim.values+weekly_std.values),
                      np.squeeze(weekly_clim.values-weekly_std.values), 
                      facecolor='steelblue', 
                      interpolate=True , 
                      alpha=.3)
    plt.ylabel(r'T ($^{\circ}C$)')
    plt.xlabel('Week of the year')
    plt.title('{} temperature'.format(siteName))
    plt.xlim([0,53])
    plt.ylim([-2,18])
    plt.grid()
    plt.legend(['1989-2018 average', year])
    fig = ax.get_figure()
    fig.set_size_inches(w=12,h=8)
    fig_name = '{}_T.png'.format(siteName)   
    fig.savefig(fig_name, dpi=200)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)
    
    plt.show()
    

# =============================================================================
# Plots the standard anomalies for a set of months. Takes in a dataframe (that has
# has been resampled and averaged by month), and a number for the lower month 
# bound and a number for the upper month bound (ie. for June - July, input 6 
# and 7). 
#
# Plots graph.
# =============================================================================
def plotAnomalies(df_monthly, lowerMonthNum, upperMonthNum, siteName):
    
    df_summer = df_monthly[(df_monthly.index.month>= lowerMonthNum) & 
                           (df_monthly.index.month<= upperMonthNum)]
    df_summer = df_summer.resample('As').mean()
    df_summer.index = df_summer.index.year
    #df_summer.to_csv('comfort_cove_thermograph_1989-2017_June-July.csv')
    
    
    ## ---- plot summer data in anomalies ---- ##
    anom = (df_summer - df_summer.mean()) / df_summer.std()
    df1 = anom[anom<0]
    df2 = anom[anom>0]
    fig = plt.figure(4)
    fig.clf()
    width = .9
    plt.bar(df1.index, np.squeeze(df1.values), width, alpha=0.8, 
                 color='steelblue')
    plt.bar(df2.index, np.squeeze(df2.values), width, bottom=0, 
                 alpha=0.8, color='indianred')
    plt.ylabel('Standardized Anomaly')
    plt.xlabel('Year')
    plt.title('{} temperature ({}-{})'.format(siteName, 
                                              cal.month_name[lowerMonthNum],
                                              cal.month_name[upperMonthNum])) 
    plt.grid()
    fig.set_size_inches(w=15,h=9)
    fig_name = '{}_anomalies.png'.format(siteName)
    #plt.annotate('data source: NCDC/NOAA', xy=(.75, .01), xycoords='figure fraction', annotation_clip=False, FontSize=12)
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)
    plt.show()
    
    
# =============================================================================
# Plots time series and regression line for a set of months. Takes in a data
# frame, the number of the lower and upper months, and the site name.  
#
# Plots graph.
# =============================================================================
def plotTSMonth(df_monthly, lowerMonthNum, upperMonthNum, siteName):
    
    #takes df_monthly and creates a new df between the given bounds
    df_series = df_monthly[(df_monthly.index.month>= lowerMonthNum) & 
                            (df_monthly.index.month<= upperMonthNum)]
    df_series = df_series.resample('As').mean()
    df_series.index = df_series.index.year


    sns.regplot(x=df_series.index , y=df_series['temperature'], 
                data = df_series)
   
    plt.xlabel("Year")
    plt.ylabel(r'T ($^{\circ}C$)')
    plt.title('{} time series ({}-{})'.format(siteName, 
                                              cal.month_name[lowerMonthNum],
                                               cal.month_name[upperMonthNum]))
    plt.grid()
    plt.show()


# =============================================================================
# plotting time series for one month out of a certain year 
# =============================================================================
def plotTSDay(df, year, month, siteName):
    df_series = df[df.index.year == year]
    
    df_day = df_series[df_series.index.month == month]
    
    fig, ax = plt.subplots()
    
    plt.plot(df_day.index, df_day['temperature'])
    
    #ax.set_xticklabels(df.index.day)
    
    plt.rc('xtick', labelsize= 5)
    plt.xlabel("Day of the month")
    plt.ylabel("Temperature")
    
    plt.title('Daily Temps for {} {} ({})'.format(cal.month_name[month],
                                                  year, siteName))
    
    plt.show()
 
        
    
    
# =============================================================================
# Plot climatology
#    
#    add default years
# =============================================================================
def plotClimatology(df, year, siteName, method = 'doy'):
    
    dfC = df
    
    if method == 'doy':
         dfC[method] = dfC.index.dayofyear
    if method == 'woy':
        dfC[method] = dfC.index.weekofyear
    
    daily_clim = dfC.groupby(method).mean() 
    daily_std = dfC.groupby(method).std() 
    df_year = dfC[dfC.index.year == year] # >
    daily_year = df_year.groupby(method).mean()
    

    #climatology time series
    ax = daily_clim.plot(linewidth=1, legend=None)
    
    #time series for inputted year
    daily_year.plot(ax=ax, linewidth=1)

    #one standard deviation away from the daily_clim value
    lowerBound = daily_clim.values - (daily_std.values)*0.5
    upperBound = daily_clim.values + (daily_std.values)*0.5
    
    plt.fill_between(daily_clim.index,
                  np.squeeze(upperBound),
                  np.squeeze(lowerBound),
                  facecolor='steelblue',
                  interpolate=True ,
                  alpha=.3)
    
    plt.grid()
    plt.rc('xtick', labelsize= 10)
    plt.rc('ytick', labelsize= 10)
    plt.xlabel("Day of the Year")
    plt.ylabel("Temperature")
    
    plt.title('{} Temps for {}'.format(year, siteName))
    ax.legend(['Climatology', year])
    
    plt.show()
    
    
    
    #first recorded day of the inputted year
    #indexDY = daily_year.index.values[0]
    
    
    # #subsets bounds so that tjey're the same length as daily_year
    # lowerBound = lowerBound[indexDY - 1:]
    # upperBound = upperBound[indexDY - 1:]
    
    
    upwellDates = []
    
    
    
    return upwellDates
    
# =============================================================================
# finding upwelling
# =============================================================================
def plotUpwells(df, year, siteName, threshold = 0.5):
    dfU = df

    dfU['doy'] = dfU.index.dayofyear

    df_year = dfU[dfU.index.year == year]
    daily_year = df_year.groupby('doy').mean()
    
    #smooths daily_year
    rolledMean = daily_year.rolling(30, center = True).mean()
    rolledStd = daily_year.rolling(30, center = True).std()

    #one std away from rolled mean 
    lowerBound = rolledMean.values - (rolledStd.values)*threshold
    upperBound = rolledMean.values + (rolledStd.values)*threshold

    #deep copy so that daily_year values arent affected
    upwellDates = copy.deepcopy(daily_year)
    
    #all values above the lower bound changed to NaN
    for i in range(len(upwellDates)):
        if((upwellDates.values[i] >= lowerBound[i]) | (np.isnan(lowerBound[i]))):
            upwellDates.values[i] = np.nan
    
    
    #plots rolled mean of daily_year
    plt.plot(rolledMean)

    
    #plots std of rolled mean of daily year
    plt.fill_between(rolledMean.index, np.squeeze(upperBound),
                  np.squeeze(lowerBound),facecolor='steelblue',
                  interpolate=True, alpha=.3)
    
    #plots daily temp average
    plt.plot(daily_year)
    
    #plots upwells
    plt.plot(upwellDates, color = 'black')
    
    plt.grid()
    plt.rc('xtick', labelsize= 10)
    plt.rc('ytick', labelsize= 10)
    plt.xlabel("Day of the Year")
    plt.ylabel("Temperature")
    
    plt.title('Upwells for {} {}'.format(siteName, year))
    #plt.legend(['Rolled Mean', year])
    
    plt.show()
    
    
    return upwellDates


# =============================================================================
# prints table with upwell info
# =============================================================================
def upwellTable(upwellDates):
    start = []
    end = []
    duration = []
    count = 0
    i = 0
    
    while i < len(upwellDates):
        if not np.isnan(upwellDates.values[i][0]):
            start.append(upwellDates.index.values[i])
            
            while(not np.isnan(upwellDates.values[i][0])):
                i += 1 
                count += 1
            
            end.append(upwellDates.index.values[i-1])
            duration.append(count)
            count = 0
            
        i += 1
       
        
    data = {'startDate': start, 'endDate': end, 'duration': duration}
    upwells = pd.DataFrame(data, 
                           columns = ['startDate', 'endDate', 'duration'])
   
    
    return upwells

# =============================================================================
# Derivative plot  
# =============================================================================
def plotDerivative(df, year, month, siteName):
    dfD = df[df.index.year == year]
    dfD = dfD[dfD.index.month == month]
    
    dfD = dfD.resample('12H').mean()
    
    dx = dfD.diff()
    dt = dfD.index.to_series().diff().dt.seconds/3600
    
    dxdt = dx['temperature']/dt
    
    dxdt.plot()
    plt.xlabel("Day of the month")
    plt.ylabel("delta temp/delta time")
    plt.title("Derivative of temp {}".format(siteName))
    plt.show()
    
       

# =============================================================================
# Takes in a dataframe, the dataframe with header info, and the site name
# and converts to a netCDF file.
#
# Returns a dataset and creates a .nc file.
# =============================================================================
def convertNetCDF(df_all, headersdf, siteName):
    
    #converting dataframe to dataset
    #dimension = datetime variable
    #data variable = temperature
    xr = xarray.Dataset.from_dataframe(df_all)
    
    #attribute values
    #creates errors when attribute is a list, so need to convert to string
    station = str(headersdf['Station'].unique())
    lat = str(headersdf['Latitude'].unique())
    long = str(headersdf['Longitude'].unique())
    startDate = str(headersdf['Start Date'].unique())
    endDate = str(headersdf['End Date'].unique())
    instDepth = str(headersdf['Inst Depth'].unique())
    fileName = str(headersdf['File Name'].unique())
    
    #sets the attributes
    xr.attrs={'Station': station,'Latitude': lat, 'Longitude': long,
              'Start Date': startDate, 'End Date': endDate, 'Inst Depth':instDepth,
              'File Name': fileName}
    
    #sets characteristics to data variable
    xr['temperature'].attrs={'units':'Celcius', 'long_name':'Temperature'}
    
    #converts to netCDF file where file name = 'siteName_netCDF.nc'
    xr.to_netcdf('{}_netCDF.nc'.format(siteName))
    
    return xr

