########################################################################################################################
##
##	GLOBAL AVERAGE FILTER
##
##	Define the filtering function:
##	Input: velocity
##	Output: index of spikes after several loops.
import numpy as np
def globalAverageFilter(U,window,data,method,writePaths_figures,VariableName):
#
##	This filter operates by estimating the mean and std of the input data and identifying
##	spikes that lie outside of the (mean +/- 2*std) range
##	This could be sensitive to poor data - perhaps revisit at a later date.
	mean = np.mean(U)
	std = np.std(U)
	spikes = (U<mean-2*std)+(U>mean+2*std)
	return spikes

#########################################################################################################################
