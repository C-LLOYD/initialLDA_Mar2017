###########################################################################################
##############		RAW DATA PROCESSING FUNCTIONS	###################################
####
####		Data is currently processed using three functions:
####
####	txtToDataFrame:			This function reads in the txt files and returns a
####					data frame consisting of probe position, sample
####					number, time stamp, residence time, and two
####					velocity components.
####
####	rawToProcessed_unweighted:	Adds additional columns onto the dataFrame such as
####					means, RMS velocities and reynolds stresses.
####
####	plotter:			...
####
###########################################################################################
####
##		Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
#
###########################################################################################
####		FUNCTION DEFINITON:	txtToDataFrame				###########
####
####		This function takes a given txt file, consisting of LDA data, and returns
####		a dataFrame for further analysis. Data frame consists of several columns:
####			column 1:	NXYZ	(file number and probe position)
####			column 2:	sampleNumber
####			column 3:	timeStamp
####			column 4:	resTime
####			column 5:	Ux
####			column 6:	Uy
####
###########################################################################################
def txtToDataFrame (fileName,writePath_dataFrames):
#	Open the input file and store the line data as 'content'
	with open(fileName) as f:
		content=f.readlines()
		content = [x.strip() for x in content]
##	Content is a list of strings, separated by \n, containing each line
##	First few lines represent header
##
##	Identify probe position:
##	Achieved by testing to see if the current line, when split by ';'
##	returns a length of 4. This length corresponds to the correct line.
#
	i=0
	while len(content[i].split(';')) is not 4:	
		i=i+1
##	Store the row containing the position and splits it by ';'					
	NXYZ = content[i].split(';')
##	Clean up the list NXYZ and store only scalars
	for i in range(len(NXYZ)):				
		NXYZ[i] = re.sub(' mm$','',NXYZ[i])	
##	Now check to find the raw data
##	Data column names are split by tabs and the first entry is always "Row#"
##	We use another while loop to find this row number and then we add one to find the data
	i=0
	for i in range(len(content)):
		if re.match(r'"Row#"',content[i].split('\t')[0]):
			index=i+1
			break
##	Now we have located this start of the data by the index variable
##	Now loop through the length of content and extract variables
##	Extract sample number, time stamp and two components of velocity
#
#	txt file arranged as:
#			column 0:	sampleNumber
#			column 1:	timeStamp
#			column 2:	resTime
#			column 3:	Ux
#			column 4:	Uy
#
	sampleNumber,timeStamp,resTime,Ux,Uy = [],[],[],[],[]
	for i in range(index,len(content)):
		sampleNumber.append(content[i].split('\t')[0])
		timeStamp.append(content[i].split('\t')[1])
		resTime.append(content[i].split('\t')[2])
		Ux.append(content[i].split('\t')[3])
		Uy.append(content[i].split('\t')[4])
#
#	Store variables as lists of floats
	sampleNumber=[float(i) for i in sampleNumber]
	timeStamp=[float(i) for i in timeStamp]
#	change units to seconds
	timeStamp = np.divide(timeStamp,1000)
#	residence time is currently in us, may need converting later
	resTime=[float(i) for i in resTime]
	Ux=[float(i) for i in Ux]
	Uy=[float(i) for i in Uy]
#
##	Change data to a dataFrame:
##	First change each list to a series, then combine the series.
##	Output the final data frame
	NXYZ = pd.Series(NXYZ)
	sampleNumber = pd.Series(sampleNumber)
	timeStamp = pd.Series(timeStamp)
	Ux = pd.Series(Ux)
	Uy = pd.Series(Uy)
	data = pd.DataFrame({'NXYZ':NXYZ,'sampleNumber':sampleNumber,'timeStamp':timeStamp,'resTime':resTime,'Ux':Ux,'Uy':Uy})
#
##	Write data frame as a 'pickle' which can be read in during plotting, if necessary.
##	File name is determined by probe position	
	data.to_pickle(writePath_dataFrames+'x_'+data.NXYZ[1]+'_z_'+data.NXYZ[3]+'_data_raw.pkl')
	return data;

###########################################################################################
####	Function Definition:		rawToProcessed_unWeighted		###########
####
####	This function calculates statistical quantities from the data frame input (as 
####	described above). The statistical quantities are added as additional columns to the
####	original data frame. This function does not include additional weighting terms to 
####	account for velocity bias.
####
####	Statistical quantities are calculated for each data entry using the previous values
####	i.e this is used to check convergence of stats.
####
####		Current stats: UxMean, UyMean, uxRMS, uyRMS, uv
####
###########################################################################################
def rawToProcessed_unweighted (data,writePath_dataFrames):
#
##	Read in data frame and convert required variables to lists: This speeds up looping
##	process.
	Ux = data.Ux.tolist()
	Uy = data.Uy.tolist()
#
##	Initialise statistics and loop through the length of the time series
	UxMean,UyMean,uxRMS,uyRMS,uv = [],[],[],[],[]
	for N in range(len(Ux)):
#		print(N, 'of', len(Ux))
#
##		initialise new values of means and calculate them based on previous N
##		append the UxMean values with updated values
		UxMeanNew, UyMeanNew = [],[]
		UxMeanNew = np.divide(sum(Ux[0:N+1]),(N+1))
		UxMean.append(UxMeanNew)
		UyMeanNew = np.divide(sum(Uy[0:N+1]),(N+1))
		UyMean.append(UyMeanNew)
#
##		calculate RMS velocities and Reynolds stresses
		uxRMS.append(np.divide(sum(np.power((Ux[0:N+1]-UxMeanNew),2)),N+1))
		uyRMS.append(np.divide(sum(np.power((Uy[0:N+1]-UyMeanNew),2)),N+1))
		uv.append(np.divide(sum((Ux[0:N+1]-UxMeanNew)*(Uy[0:N+1]-UyMeanNew)),N+1))
#
##	Print output of selected stats against global values using numpy library
##	If these quantities do not match check the code!
	print(np.mean(Ux))
	print(UxMean[-1])
	print(np.sqrt(uxRMS[-1]))
	print(np.std(Ux))
#
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	UxMean = pd.Series(UxMean)
	UyMean = pd.Series(UyMean)
	uxRMS = pd.Series(uxRMS)
	uyRMS = pd.Series(uyRMS)
	uv = pd.Series(uv)
	data['UxMean']= UxMean
	data['UyMean']= UyMean
	data['uxRMS']= uxRMS
	data['uyRMS']= uyRMS
	data['uv']= uv
####		NEEDS CHANGING : DOESN'T CONTAIN FLOW RATE!
	data.to_pickle(writePath_dataFrames+'x_'+data.NXYZ[1]+'_z_'+data.NXYZ[3]+'_data_unweighted.pkl')
	return data;

###########################################################################################
####	Function Definition:		plotter					###########
####
####	This function takes a data frame, currently just consisting of raw data, and plots
####	the convergence of the calculated statistics. This will later be adapted to plot 
####	both filtered data and un-biased velocity data.
####
####		Current stats: UxMean, UyMean, uxRMS, uyRMS, uv
####
###########################################################################################
#
##	Embedded bounds function: This is used to calculate +/- bounds for convergence and
##	used to limit axis range.
##
##	Input: 	convergedValue:	This values is the estimated converged value of the stat
##				in question. This is either taken as the final value in
##				the variable vector, or it is averaged over the last 30
##				seconds of sampling.
##		
##		scalar:		This value is the scalar to be multiplied to the converged
##				value i.e 0.01 will give (conv + 1%(conv)) as an upper 
##				bound
##
def bound(ConvergedValue,scalar):
	b = [ConvergedValue*scalar + ConvergedValue]
	return b

##
##	This script reads in two variables to be plotted, a save string for the output
##	png, a convergence method, and an axis limiter specification.
##
##	writeString:	Path of which the png is written. This also contains the first
##			part of the file name.
##
##	d:		Data frame containing probe position and time vector (required 
##			for plotting and saving)
##
##	v:		Variable for plotting (accessed via d.v)
##	
##	convMethod:	Takes input either 'MEAN' or 'END'. Conv method is selected based
##			on this.
##
##	axis:		Takes either 'FALSE' or [upper,lower]. The upper, lower limits are
##			factors of the estimated converged values.
##			
def plotter(writeString,data, time, V, convMethod, axis,xlabel,ylabel):
#	Plot the variables
	mpl.plot(time,V,color='k')
#	Set up convergence criteria
#	If None is given, provide a converged criteria but don't plot it
#	Converged is necessary for the axis limits to work
	if convMethod == None:
		converged = pd.Series.tolist(V)[-1]
	else:
#	If a MEAN converged value is wanted then average over the last 30 seconds of samples
		if convMethod == 'MEAN':
			converged = np.mean(V.loc[time[:]>270])
		else:
#	Else we simply take the last value - This is better for steady convergence behaviour
			converged = pd.Series.tolist(V)[-1]
#	ONLY PLOT THE BOUNDS IF CONVERGENCE CRITERIA IS GIVEN
##	Set up +/-5% bounds:
		plus5 = bound(converged,0.05)
		min5 = bound(converged,-0.05)
##	Set up +/-1% bounds:
		plus1 = bound(converged,0.01)
		min1 = bound(converged,-0.01)
##	Plot these additional bounds
		mpl.plot(time,plus5*len(V),color='k',linestyle='--')
		mpl.plot(time,min5*len(V),color='k',linestyle='--')
		mpl.plot(time,plus1*len(V),color='k',linestyle='-.')
		mpl.plot(time,min1*len(V),color='k',linestyle='-.')
#	Set up axis limits
	if axis == None:
		pass
	else:
		yUpper = bound(converged,axis[0])
		yLower = bound(converged,axis[1])
		mpl.axis([np.min(time),np.max(time),yLower[0],yUpper[0]])
#
	mpl.xlabel(xlabel)
	mpl.ylabel(ylabel)
	mpl.tight_layout()
##	Set up write string:
##	Take position from data file (NXYZ)
##	Take variable name from ylabel
	writePath = writeString+'x_'+data.NXYZ[1]+'_z_'+data.NXYZ[3]+'_'+ylabel+'_unweighted.png'
	mpl.savefig(writePath)
#	mpl.show()
	mpl.close()
	return

###########################################################################################