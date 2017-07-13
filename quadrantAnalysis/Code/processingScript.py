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
import matplotlib.pyplot as plt
import pandas as pd
from processingFunctions import txtToDataFrame
from FilterFunctions import Filter

#data=pd.read_pickle("../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_moving_average_weighted.pkl")
plotTimeDep = False
plotFilDep = False
writeDataFrames = False
readData = True

if writeDataFrames == True:
	fileName = "../Data/rawData/Run9_x400_fl_8_hz_300secs.000001.txt"
	data_raw = txtToDataFrame(fileName,"../Data/processedData/")
	data_fil_ma = Filter(data_raw,'movingAverageFilter','mean',10,False,"../Data/processedData/")
	print(data_fil_ma)

if readData == True:
	data = pd.read_pickle("../Data/processedData/x_400_z_1_data_filtered_moving_average.pkl")

#
##	For now, take u component of velocity
U = data.Ux.as_matrix()
tau = data.resTime.as_matrix()
u = U*tau
#N = data.sampleNumber.as_matrix()
N = len(u)
#
##	Set up bins such that they cover whole range of u
Nbins = 99
Delta = np.linspace(np.min(u),np.max(u),Nbins+1)
DeltaCentre = (Delta[0:-1]+Delta[1:])/2
n = np.zeros(len(Delta)-1)
taui = np.zeros(len(Delta)-1)
#


for i in range(len(n)-1):
	binMin = Delta[i]; binMax = Delta[i+1];
	index = map(lambda x: binMin <= x < binMax, u)
	n[i] = np.sum(index)
	if n[i] > 0:
		taui[i] = np.mean(tau[index])

binEnd = (Delta[-2],Delta[-1])
index = map(lambda x: binEnd[0] <= x < binEnd[1], u)
n[-1] = np.sum(index)
if n[i] > 0:
	taui[-1] = np.mean(tau[index])

PDF = n*taui/((Delta[1]-Delta[0])*(np.sum(taui)))
PDF = n/((Delta[1]-Delta[0])*N)
PDF = n/(N)
#plt.plot(DeltaCentre/np.max(tau),PDF)
#plt.plot(DeltaCentre,n)
#plt.hist(u/np.max(tau),bins=99,density=True)
#plt.show()
#plt.close()

#[pdf,edges] = np.histogram(U*tau/np.max(tau),bins=99,normed=True)
plt.hist(U*tau,bins=99,normed=True)
plt.show()
plt.close()





#	
#
#######################################################################################



