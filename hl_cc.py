'''


'''

import matplotlib.pyplot as plt
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
 

# plot comparing annual curve to average curve from 1989-2018
#hl_modules.plotAnnualCurve(df_all, 2007, "ComfortCove")


# plots monthly anomalies for June-July
#hl_modules.plotAnomalies(df_all, 6,7, "ComfortCove")


# time series plots for summer months 
#hl_modules.plotMonthAverage(df_all, 6, 8, "ComfortCove")
#hl_modules.plotMonthAverage(df_all, 7, 7, "ComfortCove")


#time series for everyday out of a certain month
#hl_modules.plotDailyTS(df_all, 1997, 6, "comfort cove")


#creates an error when you run both at the same time 
#hl_modules.plotClimatology(df_all, 2007, "Comfort Cove", 'woy')
#hl_modules.plotClimatology(df_all, 2007, "Comfort Cove")



#hl_modules.plotDerivative(df_all, 2007, 7, "ComfortCove")



# upwellDates = hl_modules.getUpwellDates(df_all, 2018, "Comfort Cove")
# print(hl_modules.upwellDF(upwellDates))



siteNames = ['ArnoldsCove', 'BristolsHope', 'ComfortCove', 'Hampden', 
             'Lumsden5m', 'Lumsden15m', 'Melrose5m', 'Melrose15m',
             'OldBonaventure', 'StockCove', 'UpperGullies', 'Winterton']

for site in siteNames:
    fileList = hl_modules.getFileList(site)
    df_all = hl_modules.readrpf(fileList)
    df_all = hl_modules.removeDuplicates(df_all)
    hl_modules.upwellsToCSV(df_all, site)




# netcdf
#print(hl_modules.convertNetCDF(df_all, headerInfo, "ComfortCove"))



