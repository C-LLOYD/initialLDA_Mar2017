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
from processingFunctions import simplePlotter
from processingFunctions import doublePlotter
from FilterFunctions import Filter
from processingFunctions import errorPoly
from processingFunctions import errorCompiler
##	Currently does not use loops to define file names, since there are only 6.
fileNames = [
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000001.txt",
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000002.txt",
"../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000003.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000001.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000002.txt",
"../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000003.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000001.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000002.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000003.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000004.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000005.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000006.txt",
"../Data/rawData/16hz/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000007.txt"
]
##
writePaths_dataFrames = [
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/4hz_",
"../Data/processedData/dataFrames/8hz_",
"../Data/processedData/dataFrames/8hz_",
"../Data/processedData/dataFrames/8hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_",
"../Data/processedData/dataFrames/16hz_"
]
##
writePaths_figures = [
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/4hz_",
"../Data/processedData/figures/8hz_",
"../Data/processedData/figures/8hz_",
"../Data/processedData/figures/8hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_",
"../Data/processedData/figures/16hz_"
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
"../Data/processedData/dataFrames/8hz_x_400_z_40_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-0.5_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-0.8_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-1.0_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-20.0_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-50.0_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-5.0_data_filtered_moving_average_weighted.pkl",
"../Data/processedData/dataFrames/16hz_x_400.0_z_-90.0_data_filtered_moving_average_weighted.pkl"
]


#data=pd.read_pickle("../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_moving_average_weighted.pkl")
plotTimeDep = False
plotFilDep = False
writeDataFrames = False
plotSimpleTimeDep = True
plot16hzRuns = True

#	Error testing process:
#df1 = pd.read_pickle(dataFrames_ma_filtered[0])
#df2 = pd.read_pickle(dataFrames_ma_filtered[1])
#df3 = pd.read_pickle(dataFrames_ma_filtered[2])

#errorPoly(df1,df2,df3,'string')



if plotFilDep == True:
#
##	This loop plots data for filter type dependence - the data set and variables chosen
##	are particularly 'poor'
	data_fil_raw = pd.read_pickle(dataFrames_weighted[0])
	data_fil_ga = pd.read_pickle(dataFrames_ga_filtered[0])
	data_fil_ma = pd.read_pickle(dataFrames_ma_filtered[0])
	data_fil_ps = pd.read_pickle(dataFrames_ps_filtered[0])
##
	mpl.plot(data_fil_raw.timeStamp,data_fil_raw.uyRMS,linestyle='-',color = 'k')
	mpl.plot(data_fil_ga.timeStamp,data_fil_ga.uyRMS,linestyle='-.',color = 'k')
	mpl.plot(data_fil_ma.timeStamp,data_fil_ma.uyRMS,linestyle='--',color = 'k')
	mpl.plot(data_fil_ps.timeStamp,data_fil_ps.uyRMS,linestyle=':',color = 'k')
	mpl.show()
#
	mpl.plot(data_fil_raw.timeStamp,data_fil_raw.Uy,linestyle='-',color = 'k')
	mpl.plot(data_fil_ga.timeStamp,data_fil_ga.Uy,marker='o',linestyle=' ',color = 'k')
	mpl.plot(data_fil_ma.timeStamp,data_fil_ma.Uy,marker='*',linestyle=' ',color = 'k')
	mpl.plot(data_fil_ps.timeStamp,data_fil_ps.Uy,marker='x',linestyle=' ',color = 'k')
	mpl.show()

if writeDataFrames == True:
#	for i in [6,7,8,9,10,11,12]:
#		print(fileNames[i])
#		data_raw = txtToDataFrame(fileNames[i],writePaths_dataFrames[i])
#		print('Moving Average:')
#		data_fil_ma = Filter(data_raw,'movingAverageFilter','mean',10,writePaths_figures[i],writePaths_dataFrames[i])
#		print('Phase Space:')
#		data_fil_ps = Filter(data_raw,'phaseSpaceFilter','mean',10,writePaths_figures[i],writePaths_dataFrames[i])
#		print('Global Average:')
#		data_fil_ga = Filter(data_raw,'globalAverageFilter','mean',10,writePaths_figures[i],writePaths_dataFrames[i])
#	
#		
#
#		data_unwei 		= 	rawToProcessed_unweighted(data_raw		,writePaths_dataFrames[i],'_data_unweighted.pkl')
#		data_wei 		= 	rawToProcessed_weighted(data_raw,20		,writePaths_dataFrames[i],'_data_weighted.pkl')
#		data_fil_ma_unwei 	= 	rawToProcessed_unweighted(data_fil_ma		,writePaths_dataFrames[i],'_data_filtered_moving_average_unweighted.pkl')
#		data_fil_ma_wei 	= 	rawToProcessed_weighted(data_fil_ma,20		,writePaths_dataFrames[i],'_data_filtered_moving_average_weighted.pkl')
#		data_fil_ga_unwei 	= 	rawToProcessed_unweighted(data_fil_ga		,writePaths_dataFrames[i],'_data_filtered_global_average_unweighted.pkl')
#		data_fil_ga_wei 	= 	rawToProcessed_weighted(data_fil_ga,20		,writePaths_dataFrames[i],'_data_filtered_global_average_weighted.pkl')
#		data_fil_ps_unwei 	= 	rawToProcessed_unweighted(data_fil_ps		,writePaths_dataFrames[i],'_data_filtered_phase_space_unweighted.pkl')
#		data_fil_ps_wei 	= 	rawToProcessed_weighted(data_fil_ps,20		,writePaths_dataFrames[i],'_data_filtered_phase_space_weighted.pkl')
##
#
#	Now write error dataFrames:
#	df_ma_4hz = dataFrames_ma_filtered[0:3]
#	df_ma_8hz = dataFrames_ma_filtered[3:6]
	df_ma_16hz = dataFrames_ma_filtered[6:]
	print(df_ma_16hz)
#	errorCompiler(df_ma_4hz,'../Data/processedData/dataFrames/4hz_errors.pkl')
#	errorCompiler(df_ma_8hz,'../Data/processedData/dataFrames/8hz_errors.pkl')
	errorCompiler(df_ma_16hz,'../Data/processedData/dataFrames/16hz_errors.pkl')

##	Plotting functions:
#
##	Plot mean Ux
if plotTimeDep == True:
	axis = [-0.6,0.6]
	for i in [0,1,2,3,4,5]:
		print(fileNames[i])
		data1 = pd.read_pickle(dataFrames_weighted[i])
		data2 = pd.read_pickle(dataFrames_ma_filtered[i])
		data3 = pd.read_pickle(dataFrames_ga_filtered[i])
		data4 = pd.read_pickle(dataFrames_ps_filtered[i])
#
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
if plotSimpleTimeDep == True:
#	for i in [0,1,2,3,4,5]:
#
#	print(fileNames[i])
	data_4_1 = pd.read_pickle(dataFrames_ma_filtered[0])
	data_4_15 = pd.read_pickle(dataFrames_ma_filtered[1])
	data_4_40 = pd.read_pickle(dataFrames_ma_filtered[2])
	data_8_1 = pd.read_pickle(dataFrames_ma_filtered[3])
	data_8_15 = pd.read_pickle(dataFrames_ma_filtered[4])
	data_8_40 = pd.read_pickle(dataFrames_ma_filtered[5])
#		print(data)
#
##	Plot Ux
	simplePlotter(
	v1 = data_4_1.timeStamp,	v2 = data_4_15.timeStamp,	v3 = data_4_40.timeStamp,
	v4 = data_8_1.timeStamp,	v5 = data_8_15.timeStamp,	v6 = data_8_40.timeStamp,
	u1 = data_4_1.UxMean/np.mean(data_4_1.UxMean.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.UxMean/np.mean(data_4_15.UxMean.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.UxMean/np.mean(data_4_40.UxMean.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.UxMean/np.mean(data_8_1.UxMean.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.UxMean/np.mean(data_8_15.UxMean.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.UxMean/np.mean(data_8_40.UxMean.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4 Hz, z = 1 mm',	u2Lab = '4 Hz, z = 15 mm',	u3Lab = '4 Hz, z = 40 mm',
	u4Lab = '8 Hz, z = 1 mm',	u5Lab = '8 Hz, z = 15 mm',	u6Lab = '8 Hz, z = 40 mm',
	ylabel = r'$\frac{\mu_u(t)}{\mu_{u}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_UxMeanTConvergence.png")
#
	simplePlotter(
	v1 = data_4_1.sampleNumber,	v2 = data_4_15.sampleNumber,	v3 = data_4_40.sampleNumber,
	v4 = data_8_1.sampleNumber,	v5 = data_8_15.sampleNumber,	v6 = data_8_40.sampleNumber,
	u1 = data_4_1.UxMean/np.mean(data_4_1.UxMean.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.UxMean/np.mean(data_4_15.UxMean.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.UxMean/np.mean(data_4_40.UxMean.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.UxMean/np.mean(data_8_1.UxMean.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.UxMean/np.mean(data_8_15.UxMean.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.UxMean/np.mean(data_8_40.UxMean.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4 Hz, z = 1 mm',	u2Lab = '4 Hz, z = 15mm',	u3Lab = '4Hz, z = 40mm',
	u4Lab = '8 Hz, z = 1 mm',	u5Lab = '8 Hz, z = 15mm',	u6Lab = '8Hz, z = 40mm',
	ylabel = r'$\frac{\mu_u(t)}{\mu_{u}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/4and8hz_UxMeanNConvergence.png")
#
#
##	Plot Uy
	simplePlotter(
	v1 = data_4_1.timeStamp,	v2 = data_4_15.timeStamp,	v3 = data_4_40.timeStamp,
	v4 = data_8_1.timeStamp,	v5 = data_8_15.timeStamp,	v6 = data_8_40.timeStamp,
	u1 = data_4_1.UyMean/np.mean(data_4_1.UyMean.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.UyMean/np.mean(data_4_15.UyMean.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.UyMean/np.mean(data_4_40.UyMean.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.UyMean/np.mean(data_8_1.UyMean.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.UyMean/np.mean(data_8_15.UyMean.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.UyMean/np.mean(data_8_40.UyMean.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, z = 1mm',	u2Lab = '4Hz, z = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, z = 1mm',	u5Lab = '8Hz, z = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\mu_v(t)}{\mu_{v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_UyMeanTConvergence.png")
#
	simplePlotter(
	v1 = data_4_1.sampleNumber,	v2 = data_4_15.sampleNumber,	v3 = data_4_40.sampleNumber,
	v4 = data_8_1.sampleNumber,	v5 = data_8_15.sampleNumber,	v6 = data_8_40.sampleNumber,
	u1 = data_4_1.UyMean/np.mean(data_4_1.UyMean.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.UyMean/np.mean(data_4_15.UyMean.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.UyMean/np.mean(data_4_40.UyMean.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.UyMean/np.mean(data_8_1.UyMean.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.UyMean/np.mean(data_8_15.UyMean.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.UyMean/np.mean(data_8_40.UyMean.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\mu_v(t)}{\mu_{v}(300)}$',	xlabel = r'$N$',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_UyMeanNConvergence.png")
#
#
##	Plot Uxrms
	simplePlotter(
	v1 = data_4_1.timeStamp,	v2 = data_4_15.timeStamp,	v3 = data_4_40.timeStamp,
	v4 = data_8_1.timeStamp,	v5 = data_8_15.timeStamp,	v6 = data_8_40.timeStamp,
	u1 = data_4_1.uxRMS/np.mean(data_4_1.uxRMS.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uxRMS/np.mean(data_4_15.uxRMS.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uxRMS/np.mean(data_4_40.uxRMS.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uxRMS/np.mean(data_8_1.uxRMS.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uxRMS/np.mean(data_8_15.uxRMS.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uxRMS/np.mean(data_8_40.uxRMS.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\sigma_u(t)}{\sigma_{u}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_uxRMSTConvergence.png")
#
	simplePlotter(
	v1 = data_4_1.sampleNumber,	v2 = data_4_15.sampleNumber,	v3 = data_4_40.sampleNumber,
	v4 = data_8_1.sampleNumber,	v5 = data_8_15.sampleNumber,	v6 = data_8_40.sampleNumber,
	u1 = data_4_1.uxRMS/np.mean(data_4_1.uxRMS.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uxRMS/np.mean(data_4_15.uxRMS.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uxRMS/np.mean(data_4_40.uxRMS.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uxRMS/np.mean(data_8_1.uxRMS.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uxRMS/np.mean(data_8_15.uxRMS.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uxRMS/np.mean(data_8_40.uxRMS.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\sigma_u(t)}{\sigma_{u}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/4and8hz_uxRMSNConvergence.png")
#
#
##	Plot uyRMS
	simplePlotter(
	v1 = data_4_1.timeStamp,	v2 = data_4_15.timeStamp,	v3 = data_4_40.timeStamp,
	v4 = data_8_1.timeStamp,	v5 = data_8_15.timeStamp,	v6 = data_8_40.timeStamp,
	u1 = data_4_1.uyRMS/np.mean(data_4_1.uyRMS.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uyRMS/np.mean(data_4_15.uyRMS.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uyRMS/np.mean(data_4_40.uyRMS.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uyRMS/np.mean(data_8_1.uyRMS.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uyRMS/np.mean(data_8_15.uyRMS.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uyRMS/np.mean(data_8_40.uyRMS.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\sigma_v(t)}{\sigma_{v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_uyRMSTConvergence.png")
#
	simplePlotter(
	v1 = data_4_1.sampleNumber,	v2 = data_4_15.sampleNumber,	v3 = data_4_40.sampleNumber,
	v4 = data_8_1.sampleNumber,	v5 = data_8_15.sampleNumber,	v6 = data_8_40.sampleNumber,
	u1 = data_4_1.uyRMS/np.mean(data_4_1.uyRMS.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uyRMS/np.mean(data_4_15.uyRMS.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uyRMS/np.mean(data_4_40.uyRMS.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uyRMS/np.mean(data_8_1.uyRMS.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uyRMS/np.mean(data_8_15.uyRMS.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uyRMS/np.mean(data_8_40.uyRMS.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\sigma_v(t)}{\sigma_{v}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/4and8hz_uyRMSNConvergence.png")
#
##	Plot uv
	simplePlotter(
	v1 = data_4_1.timeStamp,	v2 = data_4_15.timeStamp,	v3 = data_4_40.timeStamp,
	v4 = data_8_1.timeStamp,	v5 = data_8_15.timeStamp,	v6 = data_8_40.timeStamp,
	u1 = data_4_1.uv/np.mean(data_4_1.uv.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uv/np.mean(data_4_15.uv.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uv/np.mean(data_4_40.uv.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uv/np.mean(data_8_1.uv.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uv/np.mean(data_8_15.uv.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uv/np.mean(data_8_40.uv.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\gamma_{u,v}(t)}{\gamma_{u,v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_uvTConvergence.png")
#
	simplePlotter(
	v1 = data_4_1.sampleNumber,	v2 = data_4_15.sampleNumber,	v3 = data_4_40.sampleNumber,
	v4 = data_8_1.sampleNumber,	v5 = data_8_15.sampleNumber,	v6 = data_8_40.sampleNumber,
	u1 = data_4_1.uv/np.mean(data_4_1.uv.loc[data_4_1.timeStamp[:]>270]),
	u2 = data_4_15.uv/np.mean(data_4_15.uv.loc[data_4_15.timeStamp[:]>270]),
	u3 = data_4_40.uv/np.mean(data_4_40.uv.loc[data_4_40.timeStamp[:]>270]),
	u4 = data_8_1.uv/np.mean(data_8_1.uv.loc[data_8_1.timeStamp[:]>270]),
	u5 = data_8_15.uv/np.mean(data_8_15.uv.loc[data_8_15.timeStamp[:]>270]),
	u6 = data_8_40.uv/np.mean(data_8_40.uv.loc[data_8_40.timeStamp[:]>270]),
#	u2 = [],u5 = [],
	u1Lab = '4Hz, y = 1mm',	u2Lab = '4Hz, y = 15mm',	u3Lab = '4Hz, y = 40mm',
	u4Lab = '8Hz, y = 1mm',	u5Lab = '8Hz, y = 15mm',	u6Lab = '8Hz, y = 40mm',
	ylabel = r'$\frac{\gamma_{u,v}(t)}{\gamma_{u,v}(300)}$',	xlabel = r'$N$',	legend = True,
	writeString = "../Data/processedData/figures/4and8hz_uvNConvergence.png")
#

if plot16hzRuns == True:
#
	d0p8 = pd.read_pickle(dataFrames_ma_filtered[7])
	d1 = pd.read_pickle(dataFrames_ma_filtered[8])
	d5 = pd.read_pickle(dataFrames_ma_filtered[9])
	d20 = pd.read_pickle(dataFrames_ma_filtered[10])
	d50 = pd.read_pickle(dataFrames_ma_filtered[11])
	d90 = pd.read_pickle(dataFrames_ma_filtered[12])
#
##	Plot Ux
	simplePlotter(
	v1 = d0p8.timeStamp,	v2 = d1.timeStamp,	v3 = d5.timeStamp,
	v4 = d20.timeStamp,	v5 = d50.timeStamp,	v6 = d90.timeStamp,
	u1 = d0p8.UxMean/np.mean(d0p8.UxMean.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.UxMean/np.mean(d1.UxMean.loc[d1.timeStamp[:]>270]),
	u3 = d5.UxMean/np.mean(d5.UxMean.loc[d5.timeStamp[:]>270]),
	u4 = d20.UxMean/np.mean(d20.UxMean.loc[d20.timeStamp[:]>270]),
	u5 = d50.UxMean/np.mean(d50.UxMean.loc[d50.timeStamp[:]>270]),
	u6 = d90.UxMean/np.mean(d90.UxMean.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\mu_u(t)}{\mu_{u}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/16hz_UxMeanTConvergence.png")
#
	simplePlotter(
	v1 = d0p8.sampleNumber,	v2 = d1.sampleNumber,	v3 = d5.sampleNumber,
	v4 = d20.sampleNumber,	v5 = d50.sampleNumber,	v6 = d90.sampleNumber,
	u1 = d0p8.UxMean/np.mean(d0p8.UxMean.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.UxMean/np.mean(d1.UxMean.loc[d1.timeStamp[:]>270]),
	u3 = d5.UxMean/np.mean(d5.UxMean.loc[d5.timeStamp[:]>270]),
	u4 = d20.UxMean/np.mean(d20.UxMean.loc[d20.timeStamp[:]>270]),
	u5 = d50.UxMean/np.mean(d50.UxMean.loc[d50.timeStamp[:]>270]),
	u6 = d90.UxMean/np.mean(d90.UxMean.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\mu_u(t)}{\mu_{u}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/16hz_UxMeanNConvergence.png")
#
##	Plot Uy
	simplePlotter(
	v1 = d0p8.timeStamp,	v2 = d1.timeStamp,	v3 = d5.timeStamp,
	v4 = d20.timeStamp,	v5 = d50.timeStamp,	v6 = d90.timeStamp,
	u1 = d0p8.UyMean/np.mean(d0p8.UyMean.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.UyMean/np.mean(d1.UyMean.loc[d1.timeStamp[:]>270]),
	u3 = d5.UyMean/np.mean(d5.UyMean.loc[d5.timeStamp[:]>270]),
	u4 = d20.UyMean/np.mean(d20.UyMean.loc[d20.timeStamp[:]>270]),
	u5 = d50.UyMean/np.mean(d50.UyMean.loc[d50.timeStamp[:]>270]),
	u6 = d90.UyMean/np.mean(d90.UyMean.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\mu_v(t)}{\mu_{v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/16hz_UyMeanTConvergence.png")
#
	simplePlotter(
	v1 = d0p8.sampleNumber,	v2 = d1.sampleNumber,	v3 = d5.sampleNumber,
	v4 = d20.sampleNumber,	v5 = d50.sampleNumber,	v6 = d90.sampleNumber,
	u1 = d0p8.UyMean/np.mean(d0p8.UyMean.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.UyMean/np.mean(d1.UyMean.loc[d1.timeStamp[:]>270]),
	u3 = d5.UyMean/np.mean(d5.UyMean.loc[d5.timeStamp[:]>270]),
	u4 = d20.UyMean/np.mean(d20.UyMean.loc[d20.timeStamp[:]>270]),
	u5 = d50.UyMean/np.mean(d50.UyMean.loc[d50.timeStamp[:]>270]),
	u6 = d90.UyMean/np.mean(d90.UyMean.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\mu_v(t)}{\mu_{v}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/16hz_UyMeanNConvergence.png")
#
##	Plot uxRMS
	simplePlotter(
	v1 = d0p8.timeStamp,	v2 = d1.timeStamp,	v3 = d5.timeStamp,
	v4 = d20.timeStamp,	v5 = d50.timeStamp,	v6 = d90.timeStamp,
	u1 = d0p8.uxRMS/np.mean(d0p8.uxRMS.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uxRMS/np.mean(d1.uxRMS.loc[d1.timeStamp[:]>270]),
	u3 = d5.uxRMS/np.mean(d5.uxRMS.loc[d5.timeStamp[:]>270]),
	u4 = d20.uxRMS/np.mean(d20.uxRMS.loc[d20.timeStamp[:]>270]),
	u5 = d50.uxRMS/np.mean(d50.uxRMS.loc[d50.timeStamp[:]>270]),
	u6 = d90.uxRMS/np.mean(d90.uxRMS.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\sigma_u(t)}{\sigma_{u}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/16hz_uxRMSTConvergence.png")
#
	simplePlotter(
	v1 = d0p8.sampleNumber,	v2 = d1.sampleNumber,	v3 = d5.sampleNumber,
	v4 = d20.sampleNumber,	v5 = d50.sampleNumber,	v6 = d90.sampleNumber,
	u1 = d0p8.uxRMS/np.mean(d0p8.uxRMS.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uxRMS/np.mean(d1.uxRMS.loc[d1.timeStamp[:]>270]),
	u3 = d5.uxRMS/np.mean(d5.uxRMS.loc[d5.timeStamp[:]>270]),
	u4 = d20.uxRMS/np.mean(d20.uxRMS.loc[d20.timeStamp[:]>270]),
	u5 = d50.uxRMS/np.mean(d50.uxRMS.loc[d50.timeStamp[:]>270]),
	u6 = d90.uxRMS/np.mean(d90.uxRMS.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\sigma_u(t)}{\sigma_{u}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/16hz_uxRMSNConvergence.png")
#
##	Plot uyRMS
	simplePlotter(
	v1 = d0p8.timeStamp,	v2 = d1.timeStamp,	v3 = d5.timeStamp,
	v4 = d20.timeStamp,	v5 = d50.timeStamp,	v6 = d90.timeStamp,
	u1 = d0p8.uyRMS/np.mean(d0p8.uyRMS.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uyRMS/np.mean(d1.uyRMS.loc[d1.timeStamp[:]>270]),
	u3 = d5.uyRMS/np.mean(d5.uyRMS.loc[d5.timeStamp[:]>270]),
	u4 = d20.uyRMS/np.mean(d20.uyRMS.loc[d20.timeStamp[:]>270]),
	u5 = d50.uyRMS/np.mean(d50.uyRMS.loc[d50.timeStamp[:]>270]),
	u6 = d90.uyRMS/np.mean(d90.uyRMS.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\sigma_v(t)}{\sigma_{v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/16hz_uyRMSTConvergence.png")
#
	simplePlotter(
	v1 = d0p8.sampleNumber,	v2 = d1.sampleNumber,	v3 = d5.sampleNumber,
	v4 = d20.sampleNumber,	v5 = d50.sampleNumber,	v6 = d90.sampleNumber,
	u1 = d0p8.uyRMS/np.mean(d0p8.uyRMS.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uyRMS/np.mean(d1.uyRMS.loc[d1.timeStamp[:]>270]),
	u3 = d5.uyRMS/np.mean(d5.uyRMS.loc[d5.timeStamp[:]>270]),
	u4 = d20.uyRMS/np.mean(d20.uyRMS.loc[d20.timeStamp[:]>270]),
	u5 = d50.uyRMS/np.mean(d50.uyRMS.loc[d50.timeStamp[:]>270]),
	u6 = d90.uyRMS/np.mean(d90.uyRMS.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\sigma_v(t)}{\sigma_{v}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/16hz_uyRMSNConvergence.png")
#
##	Plot uv
	simplePlotter(
	v1 = d0p8.timeStamp,	v2 = d1.timeStamp,	v3 = d5.timeStamp,
	v4 = d20.timeStamp,	v5 = d50.timeStamp,	v6 = d90.timeStamp,
	u1 = d0p8.uv/np.mean(d0p8.uv.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uv/np.mean(d1.uv.loc[d1.timeStamp[:]>270]),
	u3 = d5.uv/np.mean(d5.uv.loc[d5.timeStamp[:]>270]),
	u4 = d20.uv/np.mean(d20.uv.loc[d20.timeStamp[:]>270]),
	u5 = d50.uv/np.mean(d50.uv.loc[d50.timeStamp[:]>270]),
	u6 = d90.uv/np.mean(d90.uv.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\gamma_{u,v}(t)}{\gamma_{u,v}(300)}$',	xlabel = r'$t$ (s)',	legend = True,
	writeString = "../Data/processedData/figures/16hz_uvTConvergence.png")
#
	simplePlotter(
	v1 = d0p8.sampleNumber,	v2 = d1.sampleNumber,	v3 = d5.sampleNumber,
	v4 = d20.sampleNumber,	v5 = d50.sampleNumber,	v6 = d90.sampleNumber,
	u1 = d0p8.uv/np.mean(d0p8.uv.loc[d0p8.timeStamp[:]>270]),
	u2 = d1.uv/np.mean(d1.uv.loc[d1.timeStamp[:]>270]),
	u3 = d5.uv/np.mean(d5.uv.loc[d5.timeStamp[:]>270]),
	u4 = d20.uv/np.mean(d20.uv.loc[d20.timeStamp[:]>270]),
	u5 = d50.uv/np.mean(d50.uv.loc[d50.timeStamp[:]>270]),
	u6 = d90.uv/np.mean(d90.uv.loc[d90.timeStamp[:]>270]),
	u1Lab = 'z = 0.8 mm',	u2Lab = 'z = 1 mm',	u3Lab = 'z = 5 mm',
	u4Lab = 'z = 20 mm',	u5Lab = 'z = 50 mm',	u6Lab = 'z = 90 mm',
	ylabel = r'$\frac{\gamma_{u,v}(t)}{\gamma_{u,v}(300)}$',	xlabel = r'$N$',	legend = False,
	writeString = "../Data/processedData/figures/16hz_uvNConvergence.png")
#



testing = False
if testing == True:
	fig = mpl.figure()
	ax = fig.add_subplot(1, 1, 1)
	# Move left y-axis and bottim x-axis to centre, passing through (0,0)
	ax.spines['left'].set_position('center')
	ax.spines['bottom'].set_position('center')
	# Eliminate upper and right axes
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	#mpl.gca().set_aspect('equal', adjustable='box')
	# Show ticks in the left and lower axes only
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	u = (d.Ux-np.mean(d.Ux))
	v = (d.Uy-np.mean(d.Uy))
	mpl.scatter(d.Ux-np.mean(d.Ux),d.Uy-np.mean(d.Uy))
	mpl.axis([-np.max(np.abs(u)),np.max(np.abs(u)),-np.max(np.abs(v)),np.max(np.abs(v))])
	mpl.show()
	


#######################################################################################



