###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the spectral analysis tests of Mar2017.
####
###########################################################################################
####
##		Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
from processingFunctions import txtToDataFrame
from FilterFunctions import Filter
#
##	Set up write/read paths for raw data and dataFrames
#
##
fileName =	"../Data/rawData/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000002.txt"
writePaths_dataFrames = "../Data/processedData/dataFrames/4Hz_"
dataFrame_filtered = "../Data/processedData/dataFrames/4Hz_x_400_z_15_data_filtered_moving_average.pkl"
##
##	First function simply filters out noisey data from raw file and saves as pandas dataframe
writeDataFrames = False
if writeDataFrames == True:
	data_raw = txtToDataFrame(fileName,writePaths_dataFrames)
	D = Filter(data_raw,'movingAverageFilter','mean',10,writePaths_dataFrames)
#
##	New function:	ReSampling of time series data
#
##	1.	Read in dataFrame
D = pd.read_pickle(dataFrame_filtered)
#
##	Test if series is in logical order
#test = min(D.timeStamp.as_matrix()[1:] - D.timeStamp.as_matrix()[0:-1])
#print(test)
#
##	If test is positive then we have a correctly ordered time series.
#
##	Need to resample the data before we calculate covariances etc.
##	We take the average time step between points as our time step.
deltaT = np.mean(D.timeStamp.as_matrix()[1:]-D.timeStamp.as_matrix()[0:-1])
t = np.arange(0,300+deltaT,deltaT)
#print(t)
#
##	Now we march through the uniform time vector, identify velocities that lie in each
##	bin, and average these. 
##	
##	We initialise by taking the 0 value as the first value in the time Series - this is required
##	in order to ensure there are no empty bins
##
Uo = D.Ux.as_matrix()
to = D.timeStamp.as_matrix()
U = np.empty(len(t))
N = np.zeros(len(t))
U[0] = Uo[0]
if t[1] < to[0]:
	N[0] = 0
else:
	N[0] = 1 
for i in range(len(t))[1:]:
#	print(t[i])
	logic = (to[:] > t[i-1]) & (to[:] < t[i])
	N[i] = sum(logic)
print(N)
	

#######################################################################################



