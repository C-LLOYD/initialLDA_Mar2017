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
from processingFunctions import doublePlotter
from phaseSpaceFilter import phaseSpaceFilter
from movingAverageFilter import movingAverageFilter
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

dataFrames_ma_filtered = [
"../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_15_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400_z_40_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_1_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_15_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400_z_40_data_filtered_moving_average_weighted.pkl"
]
#data=pd.read_pickle("../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_moving_average_weighted.pkl")
axis = [-0.6,0.6]
for i in [0,1,2,3,4,5]:
	print(fileNames[i])
#	data_raw = txtToDataFrame(fileNames[i],writePaths_dataFrames[i])
#	data_unwei = rawToProcessed_unweighted(data_raw,writePaths_dataFrames[i],'_data_unweighted.pkl')
#	data_wei = rawToProcessed_weighted(data_unwei,writePaths_dataFrames[i],'_data_weighted.pkl')
#	data_fil = movingAverageFilter(data_raw,'mean',writePaths_figures[i],writePaths_dataFrames[i])
#	data_fil_unwei = rawToProcessed_unweighted(data_fil,writePaths_dataFrames[i],'_data_filtered_moving_average_unweighted.pkl')
#	data_fil_wei = rawToProcessed_weighted(data_fil_unwei,writePaths_dataFrames[i],'_data_filtered_moving_average_weighted.pkl')
	data1 = pd.read_pickle(dataFrames_weighted[i])
	data2 = pd.read_pickle(dataFrames_ma_filtered[i])
#
##	Plotting functions:
#
##	Plot mean Ux
	doublePlotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
	N1=data1.sampleNumber,	
	N2=data1.sampleNumber,	N3=data2.sampleNumber,	N4=data2.sampleNumber,
	U1=data1.UxMean,	U2=data1.UxMean_w,	U3=data2.UxMean,	U4=data2.UxMean_w,
	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{N}$',	ylabel=r'$\mathbf{\left<U_x\right>}$ (m/s)',
	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
	writeName='MeanUx.png',legend=None)
#
#	Plot mean Uy
	doublePlotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
	N1=data1.sampleNumber,	
	N2=data1.sampleNumber,	N3=data2.sampleNumber,	N4=data2.sampleNumber,
	U1=data1.UyMean,	U2=data1.UyMean_w,	U3=data2.UyMean,	U4=data2.UyMean_w,
	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{N}$',	ylabel=r'$\mathbf{\left<U_y\right>}$ (m/s)',
	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
	writeName='MeanUy.png',legend=None)
##
###	Plot RMS Ux
	doublePlotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
	N1=data1.sampleNumber,	
	N2=data1.sampleNumber,	N3=data2.sampleNumber,	N4=data2.sampleNumber,
	U1=np.sqrt(data1.uxRMS),	U2=np.sqrt(data1.uxRMS_w),	U3=np.sqrt(data2.uxRMS),	
	U4=np.sqrt(data2.uxRMS_w),
	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{N}$',	ylabel=r'RMS($\mathbf{u_x}$) (m/s)',
	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
	writeName='RMSux.png',legend=None)
#
###	Plot RMS Uy
	doublePlotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
	N1=data1.sampleNumber,	
	N2=data1.sampleNumber,	N3=data2.sampleNumber,	N4=data2.sampleNumber,
	U1=np.sqrt(data1.uyRMS),	U2=np.sqrt(data1.uyRMS_w),	U3=np.sqrt(data2.uyRMS),	
	U4=np.sqrt(data2.uyRMS_w),
	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{N}$',	ylabel=r'RMS($\mathbf{u_y}$) (m/s)',
	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
	writeName='RMSuy.png',legend=None)
##
##	Plot Reynolds Stresses uv
	doublePlotter(writeString=writePaths_figures[i],	data=data1,	time1=data1.timeStamp,	
	time2=data1.timeStamp,	time3=data2.timeStamp,	time4=data2.timeStamp,
	N1=data1.sampleNumber,	
	N2=data1.sampleNumber,	N3=data2.sampleNumber,	N4=data2.sampleNumber,
	U1=data1.uv,	U2=data1.uv_w,	U3=data2.uv,	U4=data2.uv_w,
	convMethod='MEAN',	axis=axis,	xlabel=r'$\mathbf{N}$',
	ylabel=r'$\mathbf{ \left< u_x u_y \right> }$ (m$^2$/s$^2$)',
	U1label='Raw',	U2label='Raw, Weighted',	U3label='Filtered', U4label='Filtered, Weighted',
	writeName='uv.png',legend=None)
#


#######################################################################################



