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
##		Note : 	There are two runs, both at 15mm from the plate. One at 8Hz and one at 4Hz.
path = "../Data/rawData/4Hz/*/*.txt"
save = "../Data/processedData/4Hz_15mm_"
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
#
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.xticks(fontsize=25)
mpl.yticks(fontsize=25)
mpl.plot(data.Xposition,data.UxMean,marker = 'o',linestyle = '-',color='k',linewidth='2')
mpl.xlabel(r'$\mathbf{x}$ (mm)',fontsize=30)
mpl.ylabel(r'$\mathbf{\left<U_x\right>}$ (m/s)',fontsize=30)
mpl.savefig(save+'MeanUx.png')
mpl.close()
#
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.xticks(fontsize=25)
mpl.yticks(fontsize=25)
mpl.plot(data.Xposition,data.UyMean,marker = 'o',linestyle = '-',color='k',linewidth='2')
mpl.xlabel(r'$\mathbf{x}$ (mm)',fontsize=30)
mpl.ylabel(r'$\mathbf{\left<U_y\right>}$ (m/s)',fontsize=30)
mpl.savefig(save+'MeanUy.png')
mpl.close()
#
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.xticks(fontsize=25)
mpl.yticks(fontsize=25)
mpl.plot(data.Xposition,data.uxRMS,marker = 'o',linestyle = '-',color='k',linewidth='2')
mpl.xlabel(r'$\mathbf{x}$ (mm)',fontsize=30)
mpl.ylabel(r'RMS($\mathbf{u_x}$) (m/s)',fontsize=30)
mpl.savefig(save+'RMSux.png')
mpl.close()
#
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.xticks(fontsize=25)
mpl.yticks(fontsize=25)
mpl.plot(data.Xposition,data.uyRMS,marker = 'o',linestyle = '-',color='k',linewidth='2')
mpl.xlabel(r'$\mathbf{x}$ (mm)',fontsize=30)
mpl.ylabel(r'RMS($\mathbf{u_y}$) (m/s)',fontsize=30)
mpl.savefig(save+'RMSuy.png')
mpl.close()
#
mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.xticks(fontsize=25)
mpl.yticks(fontsize=25)
mpl.plot(data.Xposition,data.uv,marker = 'o',linestyle = '-',color='k',linewidth='2')
mpl.xlabel(r'$\mathbf{x}$ (mm)',fontsize=30)
mpl.ylabel(r'$\mathbf{u_x u_y}$ $\mathbf{(m^2/s^2)}$',fontsize=30)
mpl.savefig(save+'uv.png')
mpl.close()
#
#######################################################################################



