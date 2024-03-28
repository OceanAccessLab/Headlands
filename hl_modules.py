"""
Modules for Headlands project:
    1. Get file list from a folder
    2. Read .rpf files
    3. Extract header info
    4. Remove duplicates
    5. Plot annual curve
    6. Plot anomalies
    7. Plot monthly average with regression line
    8. Plot daily time series
    9. Plot climatology 
    10. Get upwell dates and plot upwells
    11. Get upwell dataframe
    12. Convert upwell dataframe to csv file
    13. Plot derivative curve
    14. Convert to netCDF file


Authors:
Giulia Bronzi (June-August 2020)
Frederic.Cyr@dfo-mpo.gc.ca (January/February 2022)

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
from datetime import date
from dateutil.relativedelta import relativedelta
from matplotlib.colors import from_levels_and_colors
import xarray as xr

# =============================================================================
# Takes in a folder name and creates a list of each file name.
#
# Input folders options are:
## Arnolds_Cove
## Bristols Hope
## ComfortCove
## Hampden
## Lumsden
## Melrose
## OldBonaventure
## StockCove
## UpperGullies
## Winterton
# Returns a list of file names.

# =============================================================================
  
def getFileList(folderName, pre_path='', post_path=''):
    # Check user's path
    ## if getpass.getuser() == 'cyrf0006':
    ##     #files_path = '~/data/Headlands_Trimmed/'
    ##     files_path = '~/data/Headlands/' # **Headlands is a symbolic link
    ## else:
    ##     files_path = '~/Desktop/'
    # Generate the list
    infiles = folderName + '.list'
    print(os.path.join(pre_path, folderName, post_path, '*.rpf'))
    os.system('ls ' + os.path.join(pre_path, folderName, post_path, '*.rpf') + ' > ' + infiles)
    filelist = np.genfromtxt(infiles, dtype=str)
    filelist = np.reshape(filelist, filelist.size) 
    
    return filelist
 
# =============================================================================
# Takes in the list of file names returned from getFileList() and extracts the
# data, creates an array of dataframes (one df for each file) and then 
# concatenates the array into one dataframe.
#
# Returns the dataframe 'df_all' containing all the data from a site.
# =============================================================================
    
def readrpf(filelist):     
    dfs = []
    for fname in filelist:
        try:
            df = pd.read_csv(fname, sep='\s+', header=16, names=['date','time','temperature'], encoding_errors='ignore')             
            # Parse date/time
            df['datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'])
            df.set_index('datetime', inplace=True)
            df.drop(columns=['date', 'time'], inplace=True)

            ## Deprecated::
            #df = pd.read_csv(fname, sep='\s+',  
            #                 parse_dates={'datetime': [0, 1]}, 
            #                 header=16)
            #df = df.set_index('datetime')
            #df.columns = ['temperature']
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
# Takes in the list of file names returned from getFileList() and creates a 
# dataframe containting all header information. Each row of the dataframe 
# relates to one file. 
#   
# Returns the dataframe 'headersdf'.
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
                'SiteName': siteName,
                'StartDate': startDate,
                'StartTime': startTime,
                'EndDate': endDate,
                'EndTime': endTime,
                'Latitude': latitude,
                'Longitude': longitude,
                'InstType': instType,
                'SerialNumber': serialNumber,
                'WaterDepth': waterDepth,
                'InstDepth': instDepth,
                'SamplingInterval': samplingInterval,
                'FileName': fileName }
    
    colNames = ['Station', 'SiteName', 'StartDate', 'StartTime', 'EndDate', 
                'EndTime','Latitude', 'Longitude',  'InstType', 'SerialNumber',
                'WaterDepth','InstDepth', 'SamplingInterval', 'FileName']
    
    #creates new dataframe, each column is a component of the header
    #each row represents an individual file
    headersdf = pd.DataFrame(headers, columns = colNames)
    
    #sorts by start date and end date
    headersdf = headersdf.sort_values(by=['StartDate', 'EndDate'])
    #reorders the index
    headersdf = headersdf.reset_index(drop=True)

    return headersdf


# =============================================================================
# Takes in the dataframe returned from readrpf(), finds the average temperature
# between any duplicates, removes the duplicates, and returns a dataframe. 
#
# Returns the dataframe 'df'.
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
# Takes in the dataframe returned from readrpf() and a year, and produces a 
# plot comparing the inputted year with the average temperature from 1989-2018.
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
    plt.legend(['1989-2019 average', year])
    fig = ax.get_figure()
    fig.set_size_inches(w=12,h=8)
    
    plt.show()
    

# =============================================================================
# Plots the standard anomalies for a set of months. Takes in the dataframe 
# returned from readrpf(), and a number for the lower month 
# bound and a number for the upper month bound (ie. for June - July, input 6 
# and 7). 
# 
# Plots graph.
# =============================================================================
def plotAnomalies(df_all, lowerMonthNum, upperMonthNum, siteName):
    
    # monthly average
    df_monthly = df_all.resample('M').mean()

    df_summer = df_monthly[(df_monthly.index.month>= lowerMonthNum) & 
                           (df_monthly.index.month<= upperMonthNum)]
    df_summer = df_summer.resample('As').mean()
    df_summer.index = df_summer.index.year
    
    
    #calculates anomalies using the mean and std of df_summer
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
    plt.show()
    
    
# =============================================================================
# Plots the average temperature and regression line for a set of months. Takes 
# in the data frame returned from readrpf(), the number representing the lower
# and upper months, and the site name.  
#
# Plots graph.
# =============================================================================
def plotMonthAverage(df_all, lowerMonthNum, upperMonthNum, siteName):
    
    # monthly average
    df_monthly = df_all.resample('M').mean()
    
    #takes df_monthly and creates a new df between the given bounds
    df_series = df_monthly[(df_monthly.index.month>= lowerMonthNum) & 
                            (df_monthly.index.month<= upperMonthNum)]
    df_series = df_series.resample('As').mean()
    df_series.index = df_series.index.year


    sns.regplot(x=df_series.index , y=df_series['temperature'], 
                data = df_series)
   
    plt.xlabel("Year")
    plt.ylabel(r'T ($^{\circ}C$)')
    plt.title('{} Regression Line ({}-{})'.format(siteName, 
                                              cal.month_name[lowerMonthNum],
                                               cal.month_name[upperMonthNum]))
    plt.grid()
    plt.show()


# =============================================================================
# Plots the time series for a given month and year, for the specified site.
# Takes in the data frame returned from readrpf(), a year, a month, and a site 
# name.
#
# Plots graph.
# =============================================================================
def plotDailyTS(df, year, month, siteName):
    
    df_series = df[df.index.year == year]
    
    df_day = df_series[df_series.index.month == month]
    
    fig, ax = plt.subplots()
    
    plt.plot(df_day.index, df_day['temperature'])
    
    
    plt.rc('xtick', labelsize= 8)
    plt.xlabel("Day of the month")
    plt.ylabel("Temperature")
    
    plt.title('Daily Temps for {} {} ({})'.format(cal.month_name[month],
                                                  year, siteName))
    
    plt.show()
 
# =============================================================================
# Plots the climatology curve for a given site along with a time series
# from a specific year. Takes in an cn_file, a year, a site name (for title),
# and a method ('doy' or 'woy').
#
# * Note that this is a change from original version that was taking a DataFrame
# as input
#
# The default method is 'doy', meaning the annual curve follows the day of the 
# year. The other option is 'woy', where the annual curve follows the week of
# the year.
#
# Plots graph.
#
# =============================================================================
def plotClimatology(nc_file, year, clim_years=[1991,2020], siteName=None, method = 'doy', XLIM = [0, 365]):

    # Open the netCDF file
    ds = xr.open_dataset(nc_file)  
    df = ds.temperature.to_pandas()
    
    if method == 'doy':
         df[method] = df.index.dayofyear
    elif method == 'woy':
        df[method] = df.index.weekofyear
    else:
        print('Wrong method. Please specify: doy/woy')
        return

    df_clim = df[(df.index.year>=clim_years[0]) & (df.index.year<=clim_years[1])]
    daily_clim = df_clim.groupby(method).mean() 
    daily_std = df_clim.groupby(method).std() 
    df_year = df[df.index.year == year] # >
    daily_year = df_year.groupby(method).mean()
    

    #climatology time series
    fig, ax = plt.subplots() 
    daily_clim.plot(ax=ax, linewidth=1, legend=None)
    
    #time series for inputted year
    if daily_year.size>0:
        daily_year.plot(ax=ax, linewidth=1)

    #half of one standard deviation away from the daily_clim value
    lowerBound = daily_clim.values - (daily_std.values)*0.5
    upperBound = daily_clim.values + (daily_std.values)*0.5
    
    plt.fill_between(daily_clim.index,
                  np.squeeze(upperBound),
                  np.squeeze(lowerBound),
                  facecolor='steelblue',
                  interpolate=True,
                  alpha=.3)
    
    plt.grid()
    plt.xlim(XLIM)
    plt.rc('xtick', labelsize= 10)
    plt.rc('ytick', labelsize= 10)
    plt.ylabel(r"Temperature ($\rm ^{\circ}C$)")
    if method == 'doy':
        plt.xlabel("Day of the Year")
    elif method == 'woy':    
        plt.xlabel("Week of the Year")
    plt.title('Coastal temperature for {}'.format(siteName))
    if daily_year.size>0:
        ax.legend(['Climatology', year], loc=2)
    else:
        ax.legend(['Climatology'], loc=2)

    plt.show()
    fig.set_size_inches(w=5,h=4)
    fig_name = siteName.replace(' ','').replace("'","").replace( "(","_").replace(")","") + '_' + str(year) + '.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim -bordercolor White -border 10x10 ' + fig_name + ' ' + fig_name)

    
# =============================================================================
# Takes in the dataframe returned from readrpf(), a year, a site name, and 
# a threshold (default of 0.5) for calculating the lower and upper bounds of 
# the rolled mean curve. 
#
# The default for 'plot' is true. When it is true it will produce a plot  
# showing the rolled mean temperature and the daily temperature for the 
# specified year. A temperature will be highlighted as an upwell if it is 
# below the (std*threshold) of the rolled mean. 
#
# Returns the dataframe 'upwellDates', with the dates and temperatures of 
# all upwells from the specified year and site. 
# =============================================================================
def getUpwellDates(df, year, siteName, threshold = 0.5, plot = True):

    df['doy'] = df.index.dayofyear

    df_year = df[df.index.year == year]
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

    

    if(plot == True):
        #plots rolled mean
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
        
        
    #after we make the plot we convert the doy index to an absolute date
    #resets index so we can edit 'doy' as a column
    upwellDates = upwellDates.reset_index()
    
    #converts 'doy' to yyyy-mm-dd
    upwellDates.index = pd.DatetimeIndex(upwellDates['doy'].apply(lambda x: date(year, 1, 1) + relativedelta(days=int(x)-1)))
    
    #drops the doy column
    upwellDates = upwellDates.drop(columns = ['doy'])

    #index is yyyy-mm-dd, column for temperatures
    return upwellDates


# =============================================================================
# Takes in the dataframe returned from getUpwellDates() and creates
# a dataframe with one row for each upwell. There are three columns, one for 
# start date, end date, and duration of the upwell. 
#
# Returns the dataframe 'upwells'.
# =============================================================================
def upwellDF(upwellDates):
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
# Takes in the dataframe returned from readrpf() and a site name, and finds 
# the upwell dates for every year of the specified site. It then creates a 
# .csv file that has the start date, end date and duration of every upwell.
#  
# Creates a .csv file called 'upwells_siteName'. 
# =============================================================================
def upwellsToCSV(df_all, siteName):
    
    #array of dataframes, one df for each year 
    dfArray = []
    
    for year in range(1989,2019):
        upwellDates = getUpwellDates(df_all, year, siteName, 
                                     threshold = 0.5, plot = False)
        upwell = upwellDF(upwellDates)
        dfArray.append(upwell)
        
    allUpwells = pd.concat(dfArray, axis = 0)
    
    allUpwells.to_csv(r'upwells_{}.csv'.format(siteName), index=False)
    
    
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
# Takes in the dataframe returned from readrpf(), the dataframe returned from 
# extractHeaders(), and a site name, and converts the data to a netCDF file.
#
# Returns the dataset 'xr' and creates a .nc file called siteName_netCDF.nc'.
# =============================================================================
def convertNetCDF(df_all, headersdf, site):
    
    #converting dataframe to dataset
    #dimension = datetime variable, data variable = temperature
    xr = xarray.Dataset.from_dataframe(df_all)
    
    #attribute values
    #creates errors when attribute is a list, so need to convert to string
    station = str(headersdf['Station'].unique())
    siteName = str(headersdf['SiteName'].unique())
    start = str((headersdf['StartDate'] + ' ' + headersdf['StartTime']).values)
    end = str((headersdf['EndDate'] + ' ' + headersdf['EndTime']).values)
    latitude = str(headersdf['Latitude'].unique())
    longitude = str(headersdf['Longitude'].unique())
    instType = str(headersdf['InstType'].unique())
    serialNumber = str(headersdf['SerialNumber'].unique())
    waterDepth = str(headersdf['WaterDepth'].unique())
    instDepth = str(headersdf['InstDepth'].unique())
    samplingInterval = str(headersdf['SamplingInterval'].unique())
    fileName = str(headersdf['FileName'].unique())
    
    
    #sets the attributes
    xr.attrs={'Station': station, 'SiteName': siteName, 'Start': start,
              'End': end,'Latitude': latitude,'Longitude': longitude, 
              'InstType': instType,'SerialNumber': serialNumber, 
              'WaterDepth': waterDepth,'InstDepth': instDepth,
              'SamplingInterval': samplingInterval,'FileName': fileName }
    
    
    xr['temperature'].attrs={'units':'Celcius'}
    
    
    #converts the dataset 'xr' to a netCDF file
    xr.to_netcdf('{}_netCDF.nc'.format(site))
    

# =============================================================================
# From a NetCDF file, plot monthly clim and annual average in a barplot
# Add anomaly scorecards at bottom
# usage ex:
# hl.annual_plot_anom('Arnolds_Cove.nc', "Arnold's Cove", 2021)
# =============================================================================
def annual_plot_anom(data_file, station_name, current_year, clim_years=[1991, 2020]):

    print(station_name)
    ## ---- Build the colormap ----- ##
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
    # Common parameters
    hcell, wcell = 0.5, 0.6
    hpad, wpad = 0, 0

    ## ---- Load and prepare the data ---- ##
    ds = xr.open_dataset(data_file)   
    df = ds.temperature.to_pandas()

    df_clim = df[(df.index.year>=clim_years[0]) & (df.index.year<=clim_years[1])]
    df_monthly_clim = df_clim.groupby(df_clim.index.month).mean() 
    df_monthly_std = df_clim.groupby(df_clim.index.month).std() 
    df_year = df[df.index.year == current_year] # >
    df_monthly = df_year.groupby(df_year.index.month).mean()

    df_monthly = df_monthly.squeeze()
    df_monthly_clim = df_monthly_clim.squeeze()
    df_monthly_std = df_monthly_std.squeeze()

    ## ---- plot monthly ---- ##
    ind = np.arange(len(df_monthly_clim.keys()))  # the x locations for the groups
    ind_clim = df_monthly_clim.keys()-1
    ind_year = df_monthly.keys()-1

    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind_clim - width/2, df_monthly_clim.values, width, yerr=df_monthly_std.values*.5,
                    label='1991-2020')
    if df_monthly.size>0:
        rects2 = ax.bar(ind_year + width/2, np.squeeze(df_monthly.values), width, yerr=None,
                        label=str(current_year))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(r'$\rm T(^{\circ}C$)')
    ax.set_xticks(ind)
    ax.set_xlim([-.5, 11.5])
    ax.set_ylim([-1.5, 17.5])
    ax.set_xticklabels(['J','F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax.legend(loc=2)
    ax.yaxis.grid() # horizontal lines
    plt.title(station_name)

    ## ---- Add scorecards ---- ##
    std_anom = (df_monthly - df_monthly_clim) / df_monthly_std
    colors = cmap(normal(std_anom.values))
    cell_text = [std_anom.values.round(1)]
    the_table = ax.table(cellText=cell_text,
            rowLabels=['anom.'],
            colLabels=None,
            cellColours = [colors],
            cellLoc = 'center', rowLoc = 'center',
            loc='bottom', bbox=[0.0, -0.13, 1, 0.05])
    the_table.auto_set_font_size (False)
    the_table.set_fontsize(6)

    for key, cell in the_table.get_celld().items():
        cell_text = cell.get_text().get_text() 
        if key[1] == -1:
            cell.set_linewidth(0)
            cell.set_fontsize(10)
        elif cell_text=='nan':
            cell._set_facecolor('darkgray')
            cell._text.set_color('darkgray')
            cell.set_fontsize(0)
        else:
            cell._text.set_rotation(0)

    fig.set_size_inches(w=5,h=4)
    fig_name = data_file[:-3] + '_monthly.png'
    fig.savefig(fig_name, dpi=300)
    os.system('convert -trim ' + fig_name + ' ' + fig_name)


# =============================================================================
# From a NetCDF file, return annual mean standardizes anomaly
# (e.g., useful for scorecards)
# usage ex:
# hl.annual_std_anom('Arnolds_Cove.nc')
# hl.annual_std_anom('Arnolds_Cove.nc', anom_months = [6,7,8])
# =============================================================================    
def annual_std_anom(data_file, clim_years=[1991, 2020], anom_months = [5,6,7,8,9,10]):

    ds = xr.open_dataset(data_file)   
    df = ds.temperature.to_pandas()

    # unstack (year/month) the entire time series
    df_stack = df.groupby([(df.index.year),(df.index.month)]).mean().squeeze()
    df_unstack = df_stack.unstack()

    # unstack the climatological period
    df_clim_period = df[(df.index.year>=clim_years[0]) & (df.index.year<=clim_years[1])]
    df_clim_stack = df_clim_period.groupby([(df_clim_period.index.year),
                                            (df_clim_period.index.month)]).mean().squeeze()

    # Calculate monthly clim
    df_monthly_clim = df_clim_stack.unstack().mean()
    df_monthly_std = df_clim_stack.unstack().std()

    # Calculate monthly anomalies
    monthly_anom = df_unstack - df_monthly_clim 
    monthly_stdanom = (df_unstack - df_monthly_clim) /  df_monthly_std

    # Restricts months
    monthly_stdanom = monthly_stdanom[anom_months]

    # Threshold on no. of months to consider anomaly valid (50%)
    monthly_stdanom.loc[monthly_stdanom.count(axis=1)<np.ceil(len(anom_months)/2)]=np.nan

    # annual mean std anomaly
    anom_std = monthly_stdanom[anom_months].mean(axis=1)

    # add climatological mean and std
    anom_std.at['MEAN'] = df_monthly_clim.loc[anom_months].mean()
    anom_std.at['SD'] = df_monthly_std.loc[anom_months].mean()

    return anom_std

