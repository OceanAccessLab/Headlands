"""
Modules for headlands project:
    - getting file list from a folder
    - reading .rpf files
    - removing duplicates
    - extracting header info
    - writing to netCDF

@author: giuliabronzi
"""

import os
import getpass
import pandas as pd
import numpy as np


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
# 
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


# # =============================================================================
# # Code for converting to netCDF file
# # =============================================================================
# ##converting dataframe to dataset
# xr = xarray.Dataset.from_dataframe(df_all)

# print(xr)

# #attribute values
# station = headersdf['Station'].unique()
# lat = headersdf['Latitude'].unique()
# long = headersdf['Longitude'].unique()

# #sets the attributes
# xr.attrs={'Station': station,'Latitude': lat, 'Longitude': long}

# xr['temperature'].attrs={'units':'celcius', 'long_name':'Temperature'}


# xr.to_netcdf('practice.nc')

