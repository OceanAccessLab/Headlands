'''
data in : /home/cyrf0006/data/Headlands/
process in : /home/cyrf0006/AZMP/Headlands

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
 


# monthly average
df_monthly = df_all.resample('M').mean()
# df_monthly.plot()
# plt.title("Monthly Average")
# plt.show()

#df_monthly.to_csv('comfort_cove_thermograph_1989-2017_monthly.csv')



# plot comparing annual curve to average curve from 1989-2018
#hl_modules.plotAnnualCurve(df_all, 2007, "ComfortCove")


# plots monthly anomalies for June-July
#hl_modules.plotAnomalies(df_monthly, 6,7, "ComfortCove")



# time series plots for summer months 
#hl_modules.plotTSMonth(df_monthly, 6, 8, "ComfortCove")
#hl_modules.plotTSMonth(df_monthly, 7, 7, "ComfortCove")


#time series for everyday out of a certain month
#hl_modules.plotTSDay(df_all, 1997, 6, "comfort cove")


#creates an error when you run both at the same time 
#hl_modules.findClimatology(df_all, 2007, "Comfort Cove", 'woy')
#hl_modules.plotClimatology(df_all, 2007, "Comfort Cove")


upwellDates = hl_modules.plotUpwells(df_all, 2007, "Comfort Cove")

print(hl_modules.upwellTable(upwellDates))



#hl_modules.plotDerivative(df_all, 2007, 7, "ComfortCove")



# netcdf
#print(hl_modules.convertNetCDF(df_all, headerInfo, "ComfortCove"))

