#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import pandas as pd
##
#
##	Define the filtering function:
##	Input: velocity and method
##	Output: index of spikes after several loops.
def spikeLocator(U,W):
#
##	We loop through the length of U and take neighbouring points ....
##	For the first and last few points we can't do this! - Just use a slight weighting.
##	Define averaging window:
	Umeans = np.zeros(len(U))
	Umeans[:] = np.NAN
	std = np.zeros(len(U))
	std[:] = np.NAN
#
	for i in range(W,len(U)-W+1):
		Umeans[i] = np.mean(U[i-W:i+W])
		std[i]  = np.std(U[i-W:i+W])
#
##	Test : set up range such that if Umeans-2*std < U < Umeans+2*std we use Unew = U, Unew[test] = np.nan
	spikes = (U<Umeans-2*std)+(U>Umeans+2*std)+(np.isnan(Umeans))
	return spikes

#
##	Define function
def movingAverageFilter(data,method,writePaths_figures,writePath_dataFrames):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp.as_matrix()
	s = data.sampleNumber.as_matrix()
	resT = data.resTime.as_matrix()
#
##
	UyNew = Uy
	UxNew = Ux
	Spikes = np.isnan(UxNew)
	converged = False
#	Define HALF averaging window size
	W=10
#
##	Initialise filtered velocity fields
	while ~converged:
		test = len(UxNew[~Spikes])
		XSpikes = spikeLocator(UxNew,W)
		YSpikes = spikeLocator(UyNew,W)
		Spikes = XSpikes + YSpikes + Spikes
		N_Spikes = sum(Spikes)
		print("Number of new Ux spikes: "+str(sum(XSpikes)))
		print("Number of new Uy spikes: "+str(sum(YSpikes)))
		print("Total number of spikes: "+str(N_Spikes)+" = "+
			str(N_Spikes*float(100)/len(data.Uy.as_matrix()))+"%")
		if test==len(UxNew[~Spikes]):
			break
		else:
			for i in range(W,len(UxNew)-W+1):
				if Spikes[i]==True:
					UxNew[i] = np.mean(UxNew[i-1:i+1])
					UyNew[i] = np.mean(UyNew[i-1:i+1])
#			UxNew[Spikes]=np.mean(UxNew[Spikes-1:Spikes+1])
#			UyNew[Spikes]=np.mean(UxNew[Spikes-1:Spikes+1])
#		print(UxNew)
#
	Ux = Ux[~Spikes]
	Uy = Uy[~Spikes]		
	t  = t[~Spikes]
	resT  = resT[~Spikes]
	s  = s[~Spikes]
#
#
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	Ux = pd.Series(Ux)
	Uy = pd.Series(Uy)
	t = pd.Series(t)
	s = pd.Series(s)
	resT = pd.Series(resT)
#
	data['Ux']= Ux
	data['Uy']= Uy
	data['timeStamp']= t
	data['sampleNumber']=s
	data['resTime']=resT
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
	data.to_pickle(writePath_dataFrames+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_data_filtered_moving_average.pkl')
	return data;
##
##
##
