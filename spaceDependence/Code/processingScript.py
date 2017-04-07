###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the space-dependence tests of Mar2017.
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
##				PumpSpeed,	X,	UxMean,	UyMean,	uxRMS,	uyRMS,	uv	(z=15mm)
##	5.	remove temp dataframe from storage and read in the next file
##	6.	append the new dataframe with additional rows for each txt file
##	7.	when all data is read, save the data frame and run plotting scripts
##	8.	plot all 5 variables against X, for each pump speed and Z positions (might only be one z position so check)
##
#
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
path = "../Data/rawData/8Hz/*/*.txt"
data = []
for fileName in glob.glob(path):
	tempData = txtToDataFrame(fileName)
#
##	3.	filter data
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
data = data.sort_values(by='Xposition')
#
##	Plot data
##	Note:	need to loop through both data sets, 4Hz and 8Hz
##		could plot them together for each variable
##		Also need to estimate errors due to short averaging times (30 sec)
##		Use this to plot error bars ...
mpl.plot(data.Xposition,data.UxMean)
mpl.show()
mpl.plot(data.Xposition,data.UyMean)
mpl.show()
mpl.plot(data.Xposition,data.uxRMS)
mpl.show()
mpl.plot(data.Xposition,data.uyRMS)
mpl.show()
mpl.plot(data.Xposition,data.uv)
mpl.show()
#######################################################################################



