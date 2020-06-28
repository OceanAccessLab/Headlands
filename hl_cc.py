'''
AZMP script to read Headlands thermograph data from D. Senciall's .rpf format

data in : /home/cyrf0006/data/Headlands/
process in : /home/cyrf0006/AZMP/Headlands

Frederic.Cyr@dfo-mpo.gc.ca - February 2019

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

# df_monthly.to_csv('comfort_cove_thermograph_1989-2017_monthly.csv')

hl_modules.plotTimeSeries(df_monthly, 6, 8, "ComfortCove")

#plot comparing annual curve to average curve from 1989-2018
hl_modules.plotAnnualCurve(df_all, 2018, "ComfortCove")


#plots monthly anomalies for June-July
hl_modules.plotAnomalies(df_monthly, 6,7, "ComfortCove")


#print(hl_modules.convertNetCDF(df_all, headerInfo, "ComfortCove"))







