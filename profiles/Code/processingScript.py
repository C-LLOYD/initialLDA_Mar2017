###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the profile tests of Mar2017.
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
import glob
from processingFunctions import txtToDataFrame
from processingFunctions import timeAverage
from processingFunctions import plotter
from phaseSpaceFilter import phaseSpaceFilter
#
##########################################################################################
##
##	1.	loop over all data sets ... use wild cards to identify txt files
##	2.	import txt data and save as temp dataFrame
##	3.	filter this data frame using phase space filter
##	4.	apply averaging and save a new data frame with the columns:
##				PumpSpeed,	Z,	UxMean,	UyMean,	uxRMS,	uyRMS,	uv	(X is constant for each data set - store this in name)
##	5.	remove temp dataframe from storage and read in the next file
##	6.	append the new dataframe with additional rows for each txt file
##	7.	when all data is read, save the data frame and run plotting scripts
##	8.	plot all 5 variables against X, for each pump speed and Z positions (might only be one z position so check)
##
#
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
path = "../Data/rawData/8hz/400mm/*/*.txt"
data = []
for fileName in glob.glob(path):
	tempData = txtToDataFrame(fileName)
#
##	3.	filter data
	if isinstance(tempData,pd.DataFrame):
		tempData = phaseSpaceFilter(tempData,'mean')
#
##	4.	apply averaging and append a final data series
		dataNew  = timeAverage(tempData)
		if not isinstance(data,pd.DataFrame):
			data = dataNew
		else:
			data=data.append(dataNew)
#
##	Rearrange and save dataFrame - if wanted
##	Rearrange does not work as intended : 1,10,2, ...
data = data.sort_values(by='z',ascending=1)
print(data)
#
###############	TEMP LOCATION OF LAW OF THE WALL FUNCTION - MOVE THIS LATER ################
####
####	1. 	define new variables
##	Initialise variables using the full data set and iterate to find log-law region by removing extremes of z
U = data.UxMean.as_matrix()
Y = data.z.as_matrix()
#
##
U_lam = U[0:20]
Y_lam = Y[0:20]
U_turb = U[-50:]
Y_turb = Y[-50:]
#
##	Initialise first guess at Beta
betaHat_turb = np.sum((U_turb-np.mean(U_turb))*(np.log(Y_turb)-np.mean(np.log(Y_turb))))/np.sum((np.log(Y_turb)-np.mean(np.log(Y_turb)))**2)
alphaHat_turb = np.mean(U_turb) - betaHat_turb*np.mean(np.log(Y_turb))
#
betaHat_lam = np.sum((U_lam-np.mean(U_lam))*(Y_lam-np.mean(Y_lam)))/np.sum((Y_lam-np.mean(Y_lam))**2)
alphaHat_lam = np.mean(U_lam) - betaHat_lam*np.mean(Y_lam)
##	Plot this
loglawU = np.zeros(len(U))
loglawU[:] = betaHat_turb*np.log(Y[:]) + alphaHat_turb
lamlawU = np.zeros(len(U))
lamlawU[:] = betaHat_lam*Y[:] + alphaHat_lam
##
#
mpl.semilogx(Y,U,marker = 'o',linestyle = ' ')
mpl.semilogx(Y,loglawU)
mpl.semilogx(Y,lamlawU)
mpl.axis([0,80,0,np.max(U)+0.05*np.max(U)])
mpl.show()










##	Plot data
#mpl.plot(data.z,data.UxMean,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.UyMean,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uxRMS,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uyRMS,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uv,marker='o',linestyle=' ')
#mpl.show()
########################################################################################



