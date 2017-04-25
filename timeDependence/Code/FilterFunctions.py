#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv
#
#########################################################################################################################
#
##	Define function
##	This function controls the converging process and outputs the data file at the end
##	It calls the appropriate spike locator function, dependent on whether the user chooses
##	either moving average, global average, or ellipsoid method. If ellipsoid method is used, 
##	need to specify an additional averageMethod in order to choose either MAD or MEAN methods. 
def Filter(data,filterMethod,averageMethod,window,writePaths_figures,writePath_dataFrames):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp.as_matrix()
	s = data.sampleNumber.as_matrix()
	resT = data.resTime.as_matrix()
	halfWindow = int(window/2)
#
##
	UyNew = Uy
	UxNew = Ux
	Spikes = np.isnan(UxNew)	#Spikes are initialised based on the number of NANS in UxNew - there should be zero ... 
	converged = False
##	Initialise filtered velocity fields
##
##	Define filtering method
	if filterMethod == 'movingAverageFilter':
		from movingAverageFilter import movingAverageFilter as Filter
		fileAppend = 'filtered_moving_average.pkl'
		print('Moving Average Filter ...')
	elif filterMethod == 'phaseSpaceFilter':
		from phaseSpaceFilter import phaseSpaceFilter as Filter
		fileAppend = 'filtered_phase_space.pkl'
		print('Phase Space Filter ...')
	elif filterMethod == 'globalAverageFilter':
		from globalAverageFilter import globalAverageFilter as Filter
		fileAppend = 'filtered_global_average.pkl'
		print('Global Average Filter')
	else:
		print('No valid filtering method given ...')
#
##	
	while ~converged:
		test = len(UxNew[~Spikes])
		XSpikes = Filter(UxNew,window,data,averageMethod,writePaths_figures,'Ux')
		YSpikes = Filter(UyNew,window,data,averageMethod,writePaths_figures,'Uy')
#
		Spikes = XSpikes + YSpikes + Spikes
		N_Spikes = sum(Spikes)
#		print("Number of new Ux spikes: "+str(sum(XSpikes)))
#		print("Number of new Uy spikes: "+str(sum(YSpikes)))
		print("Total number of spikes: "+str(N_Spikes)+" = "+
			str(N_Spikes*float(100)/len(data.Uy.as_matrix()))+"%")
		if test==len(UxNew[~Spikes]):
			break
#
		else:
			UxNewNew = UxNew#np.zeros(len(UxNew))
			UyNewNew = UyNew#np.zeros(len(UxNew))
			if filterMethod == 'phaseSpaceFilter':
				UxNewNew[Spikes] = np.nan
				UyNewNew[Spikes] = np.nan	
#				pass
			
			else:
				UxNewNew[Spikes] = np.nan
				UyNewNew[Spikes] = np.nan
#				for i in range(len(UxNew)):
#					if Spikes[i]==True:
#						if i==0:
#							UxNewNew[i]=UxNew[i+1]
#							UyNewNew[i]=UyNew[i+1]
#						elif i==range(len(UxNew))[-1]:
#							UxNewNew[i]=UxNew[i-1]
#							UyNewNew[i]=UyNew[i-1]
#						else:
#							UxNewNew[i] = np.mean(UxNew[i-1:i+1])
#							UyNewNew[i] = np.mean(UyNew[i-1:i+1])
#				else:
#					UxNewNew[i] = UxNew[i]
#					UyNewNew[i] = UyNew[i]
#
#		mpl.plot(UyNew,marker='o')
#		mpl.plot(UyNewNew,marker='x')
#		mpl.show()
		UxNew = UxNewNew
		UyNew = UyNewNew
#
#	Ux = Ux[~Spikes]
#	Uy = Uy[~Spikes]		
#	t  = t[~Spikes]
#	resT  = resT[~Spikes]
#	s  = s[~Spikes]
#	Spikes = None
#
#
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	Ux = pd.Series(UxNew)
	Uy = pd.Series(UyNew)
	timeStamp = pd.Series(t)
	sampleNumber = pd.Series(s)
	resTime = pd.Series(resT)
	data2 = pd.DataFrame({'NXYZ':data.NXYZ,'sampleNumber':sampleNumber,
			'timeStamp':timeStamp,'resTime':resTime,'Ux':Ux,'Uy':Uy})
#
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
	data2.to_pickle(writePath_dataFrames+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_data_'+fileAppend)
	return data2;
##
##
##
