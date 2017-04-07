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
####	plotter:			This plots two variables against each other and 
####					saves the figure to a location specified as an
####					input. Also gives the option to apply bounds to the
####					plot, such as +/- % error lines and axis limiters.
####					These both use a separate, smaller, function called
####					'bound'.
####
###########################################################################################
####
##		Initialise python
import numpy as np			#Numpy for efficient numerics
import re				#Re for matching text in strings
import matplotlib.pyplot as mpl		#matplotlib for plotting
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import pandas as pd			#Pandas for dataFrame construction
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
####		Inputs: 	String containing name and location of the LDA txt file, 
####				and a string containing location of the saved output data 
####				frame.
####
####		Outputs:	dataFrame stored in memory and saved to the input save
####				path.
####
###########################################################################################
def txtToDataFrame (fileName):
#	Open the input file and store the line data as 'content'
	with open(fileName) as f:
		content=f.readlines()
		content = [x.strip() for x in content]
##	Content is a list of strings, separated by \n, containing each line
##	First few lines represent header
##
##	Identify probe position:
##	Achieved by testing to see if the current line, when split by ';'
##	returns a length of 4.
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
##	Now loop through the length of content and extract variables:
##	sample number, time stamp, residence time, and two components of velocity
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
	data = pd.DataFrame({'NXYZ':NXYZ,'sampleNumber':sampleNumber,
			'timeStamp':timeStamp,'resTime':resTime,'Ux':Ux,'Uy':Uy})
#
	return data;

###########################################################################################
####	Function Definition:		timeAverage				###########
####
####		This function accounts for velocity bias in the data by 	
####		using the residence time as a weighting function		
####
####	This function calculates statistical quantities from the data frame input (as 
####	described above). The statistical quantities are added as additional columns to the
####	original data frame.
####
####	Statistical quantities are calculated for each data entry using the previous values
####	i.e this is used to check convergence of stats.
####
####		Inputs:		dataFrame containing raw data and a write Path for the
####				output dataFrame.
####
####		Output:		dataFrame containing raw data and addition rows for the
####				statistics below. dataFrame stored in memory and saved.
####
####		Current stats: UxMean, UyMean, uxRMS, uyRMS, uv
####
###########################################################################################
def timeAverage (dataRaw):
#
##	Read in data frame and convert required variables to lists: This speeds up looping
##	process.
	resT = dataRaw.resTime.as_matrix()
	Ux = dataRaw.Ux.as_matrix()
	Uy = dataRaw.Uy.as_matrix()
#
#
##	Initialise statistics and loop through the length of the time series.
	UxMean,UyMean,uxRMS,uyRMS,uv = [],[],[],[],[]
	UxMean = np.divide(sum(Ux*resT),sum(resT))
	UyMean = np.divide(sum(Uy*resT),sum(resT))
	uxRMS = np.sqrt(np.divide(sum(np.power((Ux-UxMean),2)*resT),sum(resT)))
	uyRMS = np.sqrt(np.divide(sum(np.power((Uy-UyMean),2)*resT),sum(resT)))
	uv = np.divide(sum((Ux-UxMean)*(Uy-UyMean)*resT),sum(resT))
##	Create new dataFrame
	Xposition = pd.Series(dataRaw.NXYZ[1])
	UxMean = pd.Series(UxMean)
	UyMean = pd.Series(UyMean)
	uxRMS = pd.Series(uxRMS)
	uyRMS = pd.Series(uyRMS)
	uv = pd.Series(uv)
#
	data = pd.DataFrame({'Xposition':Xposition,'UxMean':UxMean, 'UyMean':UyMean, 'uxRMS':uxRMS, 'uyRMS':uyRMS, 'uv':uv})
#
##	Now need to write the data
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
##				seconds of sampling (determined by user).
##		
##		scalar:		This value is the scalar to be multiplied to the converged
##				value i.e 0.01 will give (conv + 1%(conv)) as an upper 
##				bound, -0.01 will give a lower bound. 
##
def bound(ConvergedValue,scalar):
	b = [ConvergedValue*scalar + ConvergedValue]
	return b

####
####	Plotting Function:
####
##	This script reads in two variables to be plotted, a save string for the output
##	png, a convergence method, and an axis limiter specification.
##
##	writeString:	Path of which the png is written. This also contains the first
##			part of the file name.
##
##	data:		Data frame containing probe position and time vector (required 
##			for plotting and saving).
##
##	time:		Time variable, taken from data frame. 
##
##	U:		Variable for plotting. Conv. definitions are based on this.
##
##	V:		Variable for plotting.
##
##	W:		Variable for plotting.
##	
##	convMethod:	Takes input either 'MEAN' or 'END'. Conv method is selected based
##			on this.
##
##	axis:		Takes either 'FALSE' or [upper,lower]. The upper, lower limits are
##			factors of the estimated converged values.
##
##	xlabel:		x-axis label as a string.
##
##	ylabel:		y-label as a string. This is currently also used for the save name.
##	
def plotter(**kargs):
	time1   	= 	kargs['time1'];
	time2  	= 	kargs['time2'];
	time3   	= 	kargs['time3'];
	time4  	= 	kargs['time4'];	
	U1	 	= 	kargs['U1'];		
	U2 		= 	kargs['U2'];		
	U3 		= 	kargs['U3'];
	U4		= 	kargs['U4'];		
	U1label 	= 	kargs['U1label'];
	U2label	= 	kargs['U2label'];
	U3label 	= 	kargs['U3label'];		
	U4label 	= 	kargs['U4label'];
	ylabel 	= 	kargs['ylabel'];			
	xlabel 	= 	kargs['xlabel'];		
	axis 		= 	kargs['axis'];
	writeString = 	kargs['writeString'];
	convMethod 	= 	kargs['convMethod'];
	data 		= 	kargs['data'];
	writeName	=	kargs['writeName'];
	legend	=	kargs['legend'];
#	Plot the variables
	mpl.rc('text', usetex=True)
#	mpl.rc('font', family='serif')
	plot1, = mpl.plot(time1,U1,color='k',linestyle='-.',linewidth='3',label=U1label)
	if not isinstance(U2,pd.Series):
		pass
	elif not isinstance(U3,pd.Series):
		plot2, = mpl.plot(time2,U2,color='k',linestyle='-',linewidth='2',label=U2label)
		if legend == None:
			pass
		else:
			mpl.legend(handles=[plot1, plot2])
	elif not isinstance(U4,pd.Series):
		plot2, = mpl.plot(time2,U2,color='k',linestyle='-',linewidth='2',label=U2label)
		plot3, = mpl.plot(time3,U3,color='r',linestyle='-.',linewidth='3',label=U3label)
		if legend == None:
			pass
		else:
			mpl.legend(handles=[plot1, plot2, plot3])
	else:
		plot2, = mpl.plot(time2,U2,color='k',linestyle='-',linewidth='2',label=U2label)
		plot3, = mpl.plot(time3,U3,color='r',linestyle='-.',linewidth='3',label=U3label)
		plot4, = mpl.plot(time4,U4,color='r',linestyle='-',linewidth='2',label=U4label)
		if legend == None:
			pass
		else:
			mpl.legend(handles=[plot1, plot2, plot3, plot4])
#
##	Set up convergence criteria
##	If None is given, provide a converged criteria but don't plot it
##	Converged is necessary for the axis limits to work
	if convMethod == None:
		converged = pd.Series.tolist(U4)[-1]
	else:
#	If a MEAN converged value is wanted then average over the last 30 seconds of samples
		if convMethod == 'MEAN':
			converged = np.mean(U4.loc[time4[:]>270])
		else:
#	Else we simply take the last value - This is better for steady convergence behaviour
			converged = pd.Series.tolist(U4)[-1]
#	ONLY PLOT THE BOUNDS IF CONVERGENCE CRITERIA IS GIVEN
##	Set up +/-5% bounds:
		plus5 = bound(converged,0.05)
		min5 = bound(converged,-0.05)
##	Set up +/-1% bounds:
		plus1 = bound(converged,0.01)
		min1 = bound(converged,-0.01)
##	Plot these additional bounds
		mpl.plot(time4,plus5*len(U4),color='k',linestyle=':',linewidth='1.5')
		mpl.plot(time4,min5*len(U4),color='k',linestyle=':',linewidth='1.5')
		mpl.plot(time4,plus1*len(U4),color='k',linestyle='-.',linewidth='1.5')
		mpl.plot(time4,min1*len(U4),color='k',linestyle='-.',linewidth='1.5')
#	Set up axis limits
	if axis == None:
		pass
	else:
		yUpper = bound(converged,axis[0])
		yLower = bound(converged,axis[1])
		mpl.axis([np.min(time4),np.max(time4),yLower[0],yUpper[0]])
#
	mpl.rc('font', family='serif')
	mpl.xlabel(xlabel,fontsize=30)
	mpl.ylabel(ylabel,fontsize=30)
	if legend == None:
		pass
	else:
		mpl.legend(fontsize=20)
#
#	mpl.tight_layout()
	mpl.xticks(fontsize=25)
	mpl.yticks(fontsize=25)
##	Set up write string:
##	Take position from data file (NXYZ)
##	Take variable name from writeNamestr(int(float(d.NXYZ[1])))
	writePath = writeString+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_'+writeName
	print(writePath)
	mpl.savefig(writePath)
#	mpl.show()
	mpl.close()
	return

###########################################################################################
