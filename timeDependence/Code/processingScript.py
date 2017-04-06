###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the time-dependence tests of Mar2017.
####		Data is currently processed using two functions:
####
####			txtToDataFrame:	This function reads in the txt files and returns a
####					data frame consisting of probe position, sample
####					number, time stamp, residence time, and two
####					velocity components.
####
####			rawToProcessed:	Adds additional columns onto the dataFrame such as
####					means, RMS velocities and reynolds stresses.
####
####		Currently working on the plotting scripts. The modularity of the script
####		will later allow the addition of filtering operations to the raw time 
####		series.
####
###########################################################################################
####
##		Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
from processingFunctions import txtToDataFrame
from processingFunctions import rawToProcessed_unweighted
from processingFunctions import rawToProcessed_weighted
from processingFunctions import plotter
from phaseSpaceFilter import phaseSpaceFilter
##	Currently does not use loops to define file names, since there are only 6.
fileNames = [
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000001.txt",
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000002.txt",
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000003.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000001.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000002.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000003.txt"
]
##
writePaths_dataFrames = [
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/8hz_",
"../Data/processedData/dataFrames/8hz_",
"../Data/processedData/dataFrames/8hz_"
]
##
writePaths_figures = [
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/8hz_",
"../Data/processedData/figures/8hz_",
"../Data/processedData/figures/8hz_"
]
##
dataFrames_unweighted = [
"../Data/processedData/dataFrames/4hz_x_400_z_1_data_unweighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_15_data_unweighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_40_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_1_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_15_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_40_data_unweighted.pkl"
]

dataFrames_weighted = [
"../Data/processedData/dataFrames/4hz_x_400_z_1_data_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_15_data_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_40_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_1_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_15_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_40_data_weighted.pkl"
]

dataFrames_filtered = [
"../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_15_data_filtered_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_40_data_filtered_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_1_data_filtered_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_15_data_filtered_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_40_data_filtered_weighted.pkl"
]

#def autocorr(x):
#	result = np.correlate(x, x, mode='full')
#	return result[result.size/2:]

#data = pd.read_pickle(dataFrames_filtered[0])
#U = data.Uy.as_matrix()
#t = data.timeStamp.as_matrix()
#U = U[~np.isnan(U)]
#t = t[~np.isnan(t)]

#iact.append(sum(autocorr_f[autocorr_f.size/2:]/autocorr_f[autocorr_f.size/2]))
#
#u = U-np.mean(U)
#
#timeseries = (U)
#mean = np.mean(timeseries)
#timeseries -= np.mean(timeseries)
#rho = np.correlate(timeseries, timeseries, mode='full')
#rho = rho[rho.size/2:]/rho[rho.size/2]
#timeScale = np.trapz(rho,x=t)
#print(timeScale)
#rho = autocorr(u)/np.var(u)

#mpl.plot(t,rho)
#mpl.axis([0,300,0,1])
#mpl.show()

axis = [-0.6,0.6]
for i in [0,1,2,3,4,5]:
	print(fileNames[i])
	data_raw = txtToDataFrame(fileNames[i],writePaths_dataFrames[i])
#	data_unwei = rawToProcessed_unweighted(data_raw,writePaths_dataFrames[i],'_data_unweighted.pkl')
#	data_wei = rawToProcessed_weighted(data_unwei,writePaths_dataFrames[i],'_data_weighted.pkl')
	data_fil = phaseSpaceFilter(data_raw,'mean',writePaths_figures[i],writePaths_dataFrames[i])
#	data_fil_unwei = rawToProcessed_unweighted(data_fil,writePaths_dataFrames[i],'_data_filtered_unweighted.pkl')
#	data_fil_wei = rawToProcessed_weighted(data_fil_unwei,writePaths_dataFrames[i],'_data_filtered_weighted.pkl')
#	data1 = pd.read_pickle(dataFrames_weighted[i])
#	data2 = pd.read_pickle(dataFrames_filtered[i])
#
##	Plotting functions:
#
##	Plot mean Ux
#	plotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
#	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
#	U1=data1.UxMean,	U2=data1.UxMean_w,	U3=data2.UxMean,	U4=data2.UxMean_w,
#	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{t}$ (s)',	ylabel=r'$\mathbf{<U_x>}$(m/s)',
#	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
#	writeName='MeanUx.png',legend=None)
#
##	Plot mean Uy
#	plotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
#	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
#	U1=data1.UyMean,	U2=data1.UyMean_w,	U3=data2.UyMean,	U4=data2.UyMean_w,
#	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{t}$ (s)',	ylabel=r'$\mathbf{<U_y>}$ (m/s)',
#	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
#	writeName='MeanUy.png',legend=None)
#
##	Plot RMS Ux
#	plotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
#	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
#	U1=np.sqrt(data1.uxRMS),	U2=np.sqrt(data1.uxRMS_w),	U3=np.sqrt(data2.uxRMS),	
#	U4=np.sqrt(data2.uxRMS_w),
#	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{t}$ (s)',	ylabel=r'RMS($\mathbf{u_x}$) (m/s)',
#	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
#	writeName='RMSux.png',legend=None)
#
##	Plot RMS Uy
#	plotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
#	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
#	U1=np.sqrt(data1.uyRMS),	U2=np.sqrt(data1.uyRMS_w),	U3=np.sqrt(data2.uyRMS),	
#	U4=np.sqrt(data2.uyRMS_w),
#	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{t}$ (s)',	ylabel=r'RMS($\mathbf{u_y}$) (m/s)',
#	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
#	writeName='RMSuy.png',legend=None)
#
##	Plot Reynolds Stresses uv
#	plotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
#	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
#	U1=data1.uv,	U2=data1.uv_w,	U3=data2.uv,	U4=data2.uv_w,
#	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{t}$ (s)',
#	ylabel=r'$\mathbf{u_x u_y}$ $\mathbf{(m^2/s^2)}$',
#	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
#	writeName='uv.png',legend=None)




#######################################################################################



