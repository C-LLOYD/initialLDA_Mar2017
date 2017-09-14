#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
#
########################################################################################################################
##
##	MOVING AVERAGE FILTER
##
##	Define the filtering function:
##	Input: velocity and the averaging window (multiple of 2)
##	Output: index of spikes after several loops.
def movingAverageFilter(U,window,data,method,VariableName):
#	Half the window for consistency
	W = int(window/2)
#
##	We loop through the length of U and take neighbouring points ....
##	For the first and last few points we can't do this! - Just use a slight weighting.
##	Define averaging window:
	Umeans = np.zeros(len(U))
	Umeans[:] = np.NAN
	std = np.zeros(len(U))
	std[:] = np.NAN
#
##	Assume that the windows at the extremes extend to the extremes ..
	Umeans[0:W]=np.mean(U[0:2*W])
	Umeans[-W:]=np.mean(U[-2*W:])
	std[0:W]=np.std(U[0:2*W])
	std[-W:]=np.std(U[-2*W:])
	for i in range(W,len(U)-W+1):
		Umeans[i] = np.mean(U[i-W:i+W])
		std[i]  = np.std(U[i-W:i+W])
#
##	Test : set up range such that if Umeans-2*std < U < Umeans+2*std we use Unew = U, Unew[test] = np.nan
	spikes = (U<Umeans-4*std)+(U>Umeans+4*std)
	return spikes

########################################################################################################################

