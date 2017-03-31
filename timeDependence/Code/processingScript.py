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
"../Data/processedData/dataFrames/4hz_x_400.0000_z_1.0000_data_unweighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400.0000_z_15.0000_data_unweighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400.0000_z_40.0000_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_1.0000_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_15.0000_data_unweighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_40.0000_data_unweighted.pkl"
]

dataFrames_weighted = [
"../Data/processedData/dataFrames/4hz_x_400.0000_z_1.0000_data_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400.0000_z_15.0000_data_weighted.pkl",
"../Data/processedData/dataFrames/4hz_x_400.0000_z_40.0000_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_1.0000_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_15.0000_data_weighted.pkl",
"../Data/processedData/dataFrames/8hz_x_400.0000_z_40.0000_data_weighted.pkl"
]

axis = None
for i in [0,1,2,3,4,5]:
	print(fileNames[i])
	data = pd.read_pickle(dataFrames_weighted[i])
#	test(time=data.timeStamp,U=data.UxMean_w,Ulabel='Weighted')
#	data = txtToDataFrame(fileNames[i],writePaths_dataFrames[i])
#	data = rawToProcessed_unweighted(data,writePaths_dataFrames[i])
#	data = rawToProcessed_weighted(data,writePaths_dataFrames[i])
##	Plotting functions:
##	Plot mean Ux
	plotter(writeString=writePaths_figures[i],	data=data,	time=data.timeStamp,
	U=data.UxMean_w,	V=data.UxMean,	W=[],	convMethod='MEAN',	axis=axis,
	xlabel='Sampling Time',	ylabel='Mean(Ux)',	Ulabel='Weighted',	Vlabel='Raw',
	Wlabel=[])
##	Plot mean Uy
	plotter(writeString=writePaths_figures[i],	data=data,	time=data.timeStamp,
	U=data.UyMean_w,	V=data.UyMean,	W=[],	convMethod='MEAN',	axis=axis,
	xlabel='Sampling Time',	ylabel='Mean(Uy)',	Ulabel='Weighted',	Vlabel='Raw',
	Wlabel=[])
##	Plot RMS Ux
	plotter(writeString=writePaths_figures[i],	data=data,	time=data.timeStamp,
	U=data.uxRMS_w,	V=data.uxRMS,	W=[],	convMethod='MEAN',	axis=axis,
	xlabel='Sampling Time',	ylabel='RMS(Ux)',	Ulabel='Weighted',	Vlabel='Raw',
	Wlabel=[])
##	Plot RMS Uy
	plotter(writeString=writePaths_figures[i],	data=data,	time=data.timeStamp,
	U=data.uyRMS_w,	V=data.uyRMS,	W=[],	convMethod='MEAN',	axis=axis,
	xlabel='Sampling Time',	ylabel='RMS(Uy)',	Ulabel='Weighted',	Vlabel='Raw',
	Wlabel=[])
##	Plot Reynolds Stresses uv
	plotter(writeString=writePaths_figures[i],	data=data,	time=data.timeStamp,
	U=data.uv_w,	V=data.uv,	W=[],	convMethod='MEAN',	axis=axis,
	xlabel='Sampling Time',	ylabel='Reynolds Stresses, (uv)',	Ulabel='Weighted',
	Vlabel='Raw',	Wlabel=[])




#######################################################################################



