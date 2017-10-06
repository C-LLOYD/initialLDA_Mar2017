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
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as signal
from processingFunctions import txtToDataFrame
from FilterFunctions import Filter
#
##	Set up write/read paths for raw data and dataFrames
#
##
#fileName =	"../Data/rawData/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000002.txt"
fileName =	"../Data/rawData/170816_x400_fl16hz_300secs/170816_x400_fl16hz_300secs.000005.txt"

#writePaths_dataFrames = "../Data/processedData/dataFrames/4Hz_"
writePaths_dataFrames = "../Data/processedData/dataFrames/16Hz_"

#dataFrame_filtered = "../Data/processedData/dataFrames/4Hz_x_400_z_15_data_filtered_moving_average.pkl"
dataFrame_filtered = "../Data/processedData/dataFrames/16Hz_x_400_z_20_data_filtered_moving_average.pkl"

##
##	First function simply filters out noisey data from raw file and saves as pandas dataframe
writeDataFrames = True
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
#deltaT = np.mean(D.timeStamp.as_matrix()[1:]-D.timeStamp.as_matrix()[0:-1])
#t = np.arange(0,300+deltaT,deltaT)
#print(t)
#
##	Now we march through the uniform time vector, identify velocities that lie in each
##	bin, and average these. 
##	
##	We initialise by taking the 0 value as the first value in the time Series - this is required
##	in order to ensure there are no empty bins
##
U = D.Ux.as_matrix()#[:10]
w = D.resTime.as_matrix()#[:10]
u = (U-np.mean(U))
t = D.timeStamp.as_matrix()#[:10]

#f = np.linspace(0.01, 100, 10000)
#pgram = signal.lombscargle(t, u, f)
#plt.subplot(2, 1, 1)
#plt.plot(t, U, 'b+')
#plt.subplot(2, 1, 2)
#plt.loglog(f, pgram)#np.sqrt(4*(pgram/normval)))
#plt.show()


#s = t[1:] - t[0:-1]
#dTau = np.mean(t[1:]-t[0:-1])
#tau = np.arange(0, np.max(t),dTau)
#print(tau)
#print(t)
#Rav = np.empty(len(tau))
#print(Rav)
#R = np.zeros([len(tau),len(u)])
#for k in range(len(Rav)):
#	print(int(k/len(Rav)*100),'%')#
#	i = []
#	for i in range(len(u)):
#		index = []
#		j = i
#		entry = np.abs((t[j] - t[i])/dTau - k) < 0.5
#		if entry == False:
#			if j+1 < len(u):
#				j = j + 1
#				entry = np.abs((t[j] - t[i])/dTau - k) < 0.5
##
#		print(entry)
#		while entry == True:
#			index = np.append(index,int(j))
##			print(index)
##			print(j)
#			if j+1 < len(u):
#				j = j + 1
#				entry = np.abs((t[j] - t[i])/dTau - k) < 0.5
#			else:
#				break
##			index[:] = int(index[:])
##		print(entry)
#
#		uj = []
#		bj = []
#		l = []
#		for l in range(len(index)):
#			uj = np.append(uj,u[int(index[l])]*u[i])
#			print(uj)
###			If we include fuzzy method then add it in here
#			bj = np.append(bj,
#			print(k,i,l)
	#	temp = uj*u[i]
#		R[k,i] = np.mean(uj)
#		print(R)
#	Rav[k] = np.mean(R[k][:])
#from testFunc import *

#w,atw,Tc,Fs,meanrm,locnor=False

#tau,R,f,S = fuzzy(t,U,w,False,500,1000,meanrm=True,fuzzy=True)
#plt.plot(tau,R)
#plt.show()

#plt.loglog(f,S)
#plt.show()




#		for j in range(len(u))[i:]:
#
##		In this loop we just want to flag up which entries are required for the b matrix
#			entry = np.abs((t[j] - t[i])/dTau - k) < 0.5
#			if entry == False:
#				break




#		print(index)
	#	index[t<t[i]] = False
#		print(i)
#		R[k][i] = np.sum(u[i]*u[index])/(np.sum(index)+1e-30)
#		print(R[k][i])		
#
#
#		for j in range(i,len(u)):
#			print(j)
#			index = np.abs((t[j] - t[i])/dTau - k) < 0.5
#			print(index)
#			tempUj = u[index]
#			R[k][i] = np.sum(u[i]*u[index])/np.sum(index)
#	Rav[k] = np.mean(R[k][:])
#	print(Rav)



#tempE = 2*np.fft.fft(Rav)
#tempf = np.fft.fftfreq(tau.shape[-1],d=dTau)
#index = tempf >= 0
#E = tempE[index]
#f = tempf[index]
#plt.loglog(f, E.real)
#plt.show()

##########################
##	RESAMPLING TESTS
dTau = np.mean(t[1:]-t[0:-1])
tau = np.arange(0, np.max(t),dTau)
print(dTau,1/dTau)
#
##	Define regular sample velocity
Us = np.zeros(len(tau))
Us[:] = np.nan
##
##	We must loop through length of tau and identify which velocities to average over
previousIndex = t[:] == t[0]
for i in range(len(tau)):
	index = []
#	print(previousIndex)
#
##	1.	identify U values that lie in the range (t < abs(tau[i])
	test = np.abs(tau[i]-t[:]) < dTau/2
	if np.sum(test) == 0:
		index[:] = previousIndex[:]
	else:
		index[:] = test[:]
		previousIndex[:] = test[:]
	Us[i] = np.mean(U[index])
#
##	This estimate is zeroth order accurate .... 
##	We could use a weighting method and pick the closest 5 points ...?
 

#plt.plot(tau,np.zeros(len(tau)),'x')
#plt.plot(t,np.zeros(len(t)),'o')
#plt.plot(tau,Us)
#plt.scatter(t,U)

#	lTau = tau[i]+dTau
#	lt = np.abs(t[:]
	

def autocorr(x):
	result = np.correlate(x, x, mode='full')
	result2 = result[:-1]
	return result[int(result.size/2):]

Rs  = autocorr(Us-np.mean(Us))
tempEs = 2*abs(np.fft.fft(Rs))
tempfs = np.fft.fftfreq(tau.shape[-1],d=dTau)
index = tempfs >= float(1e-5)
Es = tempEs[index]
fs = tempfs[index]
#
##	Add in a new function that smooths out the ACF

ndot = 1/dTau
c = (np.exp(-ndot))/(1-np.exp(-ndot))**2
RsTrue = np.zeros(len(Rs)-1)
RsTrue[0] = Rs[0]
for i in range(len(RsTrue)):
	RsTrue[i] = (2*c+1)*Rs[i] - c*(Rs[i-1] + Rs[i+1])

tempfsTrue =  np.fft.fftfreq(tau[:-1].shape[-1],d=dTau)
tempEsTrue =  2*abs(np.fft.fft(RsTrue))
index = tempfsTrue >= float(1e-5)
EsTrue = tempEsTrue[index]
fsTrue = tempfsTrue[index]
#EsTrue = Es*(1+fs**2/(1/(dTau**2)))-2*np.std(u)**2/((1/(dTau**3))*)#??????????? TAylor microscale??

plt.subplot(2, 3, 1)
plt.plot(tau, Rs, 'b+')
plt.xlabel('time,tau')
plt.ylabel('Autocorrelation, Rs')
plt.subplot(2, 3, 4)
plt.loglog(fs, Es)#np.sqrt(4*(pgram/normval)))
plt.xlabel('Frequency,fs')
plt.ylabel('Energy, Es')
####
plt.subplot(2, 3, 2)
plt.plot(tau[:-1], RsTrue, 'b+')
plt.xlabel('time,tau')
plt.ylabel('Autocorrelation, RsTrue')
plt.subplot(2, 3, 5)
plt.loglog(fsTrue, EsTrue)#np.sqrt(4*(pgram/normval)))
plt.xlabel('Frequency,fsTrue')
plt.ylabel('Energy, EsTrue')
#plt.show()
R = autocorr(u)
normval = t.shape[0]	
#f = np.linspace(0.01, 30, 1000)
E = 2*abs(signal.lombscargle(t, R, fs))
plt.subplot(2, 3, 3)
plt.plot(t, R, 'b+')
plt.xlabel('time,t')
plt.ylabel('Autocorrelation, R')
plt.subplot(2, 3, 6)
plt.loglog(fs, E)#np.sqrt(4*(pgram/normval)))
plt.xlabel('Frequency,f')
plt.ylabel('Energy, E')
plt.show()

#######################################################################################



