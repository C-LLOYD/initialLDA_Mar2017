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
from mpl_toolkits.axes_grid1 import host_subplot
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
	for i in [0,1,2,3,4,5,6,7,8,9,10]:#range(len(content)):
#		print((content[i].split()[0])=="\"AT",content[i].split()[0])
		if re.match(r'"AT',content[i].split()[0]) or re.match(r'"Row#"',content[i].split('\t')[0]):
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
		if len(content[i].split('\t')) == 5:
			sampleNumber.append(content[i].split('\t')[0])
			timeStamp.append(content[i].split('\t')[1])
			resTime.append(content[i].split('\t')[2])
			Ux.append(content[i].split('\t')[3])
			Uy.append(content[i].split('\t')[4])
		else:
			if len(sampleNumber) == 0:
				sampleNumber = [1]
			else:
				sampleNumber.append(sampleNumber[-1]+1)
			timeStamp.append(content[i].split('\t')[0])
			resTime.append(content[i].split('\t')[1])
			Ux.append(content[i].split('\t')[2])
			Uy.append(content[i].split('\t')[3])
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
##	Write data frame as a 'pickle' which can be read in during plotting, if necessary.
##	File name is determined by probe position	
	data.to_pickle(writePath_dataFrames+'x_'+str(abs(float(data.NXYZ[1])))+'_z_'+str(float(data.NXYZ[3]))+'_data_raw.pkl')
	return data;

###########################################################################################
####	Function Definition:		rawToProcessed_unweighted		###########
####
####	This function calculates statistical quantities from the data frame input (as 
####	described above). The statistical quantities are added as additional columns to the
####	original data frame. This function does not include additional weighting terms to 
####	account for velocity bias.
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
def rawToProcessed_unweighted (data,writePath_dataFrames,fileAppend):
#
##	Read in data frame and convert required variables to lists: This speeds up looping
##	process.
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
#
##	Initialise statistics and loop through the length of the time series.
	UxMean,UyMean,uxRMS,uyRMS,uv = [],[],[],[],[]
	for N in range(len(Ux)):
#		Print the data entry number if monitoring progress.
#		print(N, 'of', len(Ux))
#
##		initialise new values of means and calculate them based on previous N.
##		append the UxMean values with updated values.
		UxMeanNew, UyMeanNew = [],[]
		UxMeanNew = np.divide(sum(Ux[0:N+1]),(N+1))
		UxMean.append(UxMeanNew)
		UyMeanNew = np.divide(sum(Uy[0:N+1]),(N+1))
		UyMean.append(UyMeanNew)
#
##		calculate RMS velocities and Reynolds stresses.
		uxRMS.append(np.divide(sum(np.power((Ux[0:N+1]-UxMeanNew),2)),N+1))
		uyRMS.append(np.divide(sum(np.power((Uy[0:N+1]-UyMeanNew),2)),N+1))
		uv.append(np.divide(sum((Ux[0:N+1]-UxMeanNew)*(Uy[0:N+1]-UyMeanNew)),N+1))
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
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
	data.to_pickle(writePath_dataFrames+'x_'+str(float(data.NXYZ[1]))+'_z_'+str(abs(float(data.NXYZ[3])))+fileAppend)
	return data;


###########################################################################################
####	Function Definition:		rawToProcessed_weighted			###########
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
def rawToProcessed_weighted (data,averagingTime,writePath_dataFrames,fileAppend):
#
##	Read in data frame and convert required variables to lists: This speeds up looping
##	process.
	resT = data.resTime.as_matrix()
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp
#	t = np.array(data.timeStamp.as_matrix())
#
#
##	Initialise statistics and loop through the length of the time series.
	UxMean_w,UyMean_w,uxRMS_w,uyRMS_w,uv_w = [],[],[],[],[]
	for N in range(len(Ux)):
#		Print the data entry number if monitoring progress.
#		print(N, 'of', len(Ux))
#
##		initialise new values of means and calculate them based on previous N.
##		append the UxMean values with updated values.
		UxMeanNew, UyMeanNew = [],[]
		UxMeanNew = np.divide(sum(Ux[0:N+1]*resT[0:N+1]),sum(resT[0:N+1]))
		UxMean_w.append(UxMeanNew)
		UyMeanNew = np.divide(sum(Uy[0:N+1]*resT[0:N+1]),sum(resT[0:N+1]))
		UyMean_w.append(UyMeanNew)
#
##		calculate RMS velocities and Reynolds stresses.
		uxRMS_w.append(np.sqrt(np.divide(sum(np.power((Ux[0:N+1]-UxMeanNew),2)*resT[0:N+1]),sum(resT[0:N+1]))))
		uyRMS_w.append(np.sqrt(np.divide(sum(np.power((Uy[0:N+1]-UyMeanNew),2)*resT[0:N+1]),sum(resT[0:N+1]))))
		uv_w.append(np.divide(sum((Ux[0:N+1]-UxMeanNew)*(Uy[0:N+1]-UyMeanNew)*resT[0:N+1]),sum(resT[0:N+1])))
#
##	Calculate error associated with only averaging for the averagingTime specified
##	Capture indicies
#
	error_UxMean_w	=	[]
	error_UyMean_w	=	[]
	error_uxRMS_w	=	[]
	error_uyRMS_w	=	[]
	error_uv_w	=	[]
#
	for i in range(int(300-averagingTime)):
		UxNew = data.Ux.loc[(i < t[:]) & (t[:] < averagingTime + i)].as_matrix()
		UyNew = data.Uy.loc[(i < t[:]) & (t[:] < averagingTime + i)].as_matrix()
		resTime = data.resTime.loc[(i < t[:]) & (t[:] < averagingTime + i)].as_matrix()
		mean_UxNew = 	np.mean(UxNew)
		mean_UyNew = 	np.mean(UyNew)
		RMS_ux	=	np.sqrt(np.divide(sum(np.power((UxNew-mean_UxNew),2)*resTime),sum(resTime)))
		RMS_uy	=	np.sqrt(np.divide(sum(np.power((UyNew-mean_UyNew),2)*resTime),sum(resTime)))
		uv 	=	np.divide(sum((UxNew-mean_UxNew)*(UyNew-mean_UyNew)*resTime),sum(resTime))
#
		error_UxMean_w.append(np.abs(mean_UxNew-UxMean_w[-1]))
		error_UyMean_w.append(np.abs(mean_UyNew-UyMean_w[-1]))
		error_uxRMS_w.append(np.abs(RMS_ux-uxRMS_w[-1]))
		error_uyRMS_w.append(np.abs(RMS_uy-uyRMS_w[-1]))
		error_uv_w.append(np.abs(uv-uv_w[-1]))
#
	print(np.mean(error_UxMean_w)*100/UxMean_w[-1],
			np.mean(error_UyMean_w)*100/UyMean_w[-1],
			np.mean(error_uxRMS_w)*100/uxRMS_w[-1],
			np.mean(error_uyRMS_w)*100/uyRMS_w[-1],
			np.mean(error_uv_w)*100/uv_w[-1])
#	mpl.plot((error_UxMean_w/UxMean_w[-1]),linestyle=' ',marker='o')
#	mpl.show()
	error_UxMean_w	=	np.mean(error_UxMean_w)
	error_UyMean_w	=	np.mean(error_UyMean_w)
	error_uxRMS_w	=	np.mean(error_uxRMS_w)
	error_uyRMS_w	=	np.mean(error_uyRMS_w)
	error_uv_w	=	np.mean(error_uv_w)	
#
#
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	UxMean_w = pd.Series(UxMean_w)
	UyMean_w = pd.Series(UyMean_w)
	uxRMS_w = pd.Series(uxRMS_w)
	uyRMS_w = pd.Series(uyRMS_w)
	uv_w = pd.Series(uv_w)
	data['UxMean']= UxMean_w
	data['UyMean']= UyMean_w
	data['uxRMS']= uxRMS_w
	data['uyRMS']= uyRMS_w
	data['uv']= uv_w
#
##	Current code adds the errors to each individual dataframe for each XYZ location.
##	Separate code will combine these and create a compiled error dataframe.
	data['error_UxMean'] = error_UxMean_w
	data['error_UyMean'] = error_UyMean_w
	data['error_uxRMS'] = error_uxRMS_w
	data['error_uyRMS'] = error_uyRMS_w
	data['error_uv'] = error_uv_w
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
	data.to_pickle(writePath_dataFrames+'x_'+str(float(data.NXYZ[1]))+'_z_'+str(abs(float(data.NXYZ[3])))+fileAppend)
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
			mpl.legend(handles=[plot1, plot2, plot3, plot4]),
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
	writePath = writeString+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(abs(float(data.NXYZ[3]))))+'_'+writeName
	print(writePath)
	mpl.savefig(writePath)
#	mpl.show()
	mpl.close()
	return

###########################################################################################
####	Function Definition:		Double plotter				###########
####
####		Temporarily two plotters : This one plots two X axis for time and sampleNumber
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
def doublePlotter(**kargs):
	time1   	= 	kargs['time1'];
	time2  	= 	kargs['time2'];
	time3   	= 	kargs['time3'];
	time4  	= 	kargs['time4'];
	N1  	 	= 	kargs['N1'];
	N2  		= 	kargs['N2'];
	N3   		= 	kargs['N3'];
	N4  		= 	kargs['N4'];	
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
#
#	Plot the variables
#	mpl.rc('text', usetex=True)
	mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	mpl.rc('text', usetex=True)
#	mpl.rc('text.usetex':True ,'font.family':'serif', 'font.size':'15')
#
	fig = mpl.figure()
	ax1 = fig.add_subplot(1,1,1)
	ax1.plot(time1,U1,color='k',linestyle='-.',linewidth='3')
	ax1.plot(time2,U2,color='k',linestyle='-',linewidth='2')
	ax1.plot(time3,U3,color='r',linestyle='-.',linewidth='3')
	ax1.plot(time4,U4,color='r',linestyle='-',linewidth='2')
	ax1.set_xlabel('t (s)',fontsize='30')
	ax1.set_ylabel(ylabel,fontsize='30')
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
		ax1.plot(time4,plus5*len(U4),color='k',linestyle=':',linewidth='1.5')
		ax1.plot(time4,min5*len(U4),color='k',linestyle=':',linewidth='1.5')
		ax1.plot(time4,plus1*len(U4),color='k',linestyle='-.',linewidth='1.5')
		ax1.plot(time4,min1*len(U4),color='k',linestyle='-.',linewidth='1.5')
#	Set up axis limits
	if axis == None:
		pass
	else:
		yUpper = bound(converged,axis[1])
		yLower = bound(converged,axis[0])
		if yUpper > yLower:
			ax1.axis([np.min(time4),np.max(time4),yLower[0],yUpper[0]])
		else:
			ax1.axis([np.min(time4),np.max(time4),yUpper[0],yLower[0]])
#
	if legend == None:
		pass
	else:
		ax1.legend()
#
#	mpl.tight_layout()

#	mpl.axes().xaxis.grid(True,linestyle='-')
	ax2 = ax1.twiny()
	ax2.set_xlabel('N',fontsize='30')
	N = len(N1)
	ax2.set_xticks([1,2,3,4,5,6])
	ax2.set_xticklabels([int(N/6),int(2*N/6),int(3*N/6),int(4*N/6),int(5*N/6),int(N)])
#	mpl.xticks(fontsize=25)
#	mpl.yticks(fontsize=25)
##	Set up write string:
##	Take position from data file (NXYZ)
##	Take variable name from writeNamestr(int(float(d.NXYZ[1])))
	writePath = writeString+'x_'+str(float(data.NXYZ[1]))+'_z_'+str(float(data.NXYZ[3]))+'_'+writeName
	print(writePath)
	mpl.savefig(writePath)
#	mpl.show()
	mpl.close()
	return

###########################################################################################
##
##		Error polynomial
##
##		this function reads in three dataFrames (for three z positions) and a write name and creates an
##		output data frame consisting of 3 coefficients for a polynomial fit.
##		This is done for each variable: UxMean, UyMean, uxRMS, uyRMS and uv.
##
###########################################################################################
#
##	Load modules
import numpy as np
import pandas as pd
#
##	Define function: construct a quadratic out of three points,z, andcorresponding errors,E
##	This function does work but is not appropriate. Keep it incase required at a later date
def coeffCalc(z,E):
	LHS = E[0] - E[1] + (z[0]-z[1])*(E[1]-E[2])/(z[2]-z[1])
	RHS = z[0]**2 - z[1]**2 - (z[0]-z[1])*(z[2]**2 - z[1]**2)/(z[2]-z[1])
	C0 = LHS/RHS
	C1 = (E[0]-E[1]-C0*(z[0]**2-z[1]**2))/(z[0]-z[1])
	C2 = E[1]-C0*z[1]**2 - C1*z[1]
#
	print(E[0],C0*z[0]**2+C1*z[0]+C2)
	print(E[1],C0*z[1]**2+C1*z[1]+C2)
	print(E[2],C0*z[2]**2+C1*z[2]+C2)
#
	C = [C0,C1,C2]
	return C
#
##	Function:	Currently not useful: It reads in error data and fits a quadratic BUT
##			really we need more data points. Future work will simply use linear
##			interpolation.
def errorPoly(df1,df2,df3,writeName):
	z = [float(df1.NXYZ[3]),float(df2.NXYZ[3]),float(df3.NXYZ[3])]
	error_UxMean = [df1.error_UxMean[0]/np.mean(df1.UxMean),df2.error_UxMean[0]/np.mean(df2.UxMean),df3.error_UxMean[0]/np.mean(df3.UxMean)]
	error_UyMean = [df1.error_UyMean[0]/np.mean(df1.UyMean),df2.error_UyMean[0]/np.mean(df2.UyMean),df3.error_UyMean[0]/np.mean(df3.UyMean)]
	error_uxRMS = [df1.error_uxRMS[0],df2.error_uxRMS[0],df3.error_uxRMS[0]]
	error_uyRMS = [df1.error_uyRMS[0],df2.error_uyRMS[0],df3.error_uyRMS[0]]
	error_uv = [df1.error_uv[0],df2.error_uv[0],df3.error_uv[0]]
#
##	Now use coeffCalc to calculate coeffs for each error variable
	C_UxMean = coeffCalc(z,error_UxMean)
	C_UyMean = coeffCalc(z,error_UyMean)
	C_uxRMS = coeffCalc(z,error_uxRMS)
	C_uyRMS = coeffCalc(z,error_uyRMS)
	C_uv = coeffCalc(z,error_uv)
#
##
	X = np.linspace(0,100,100)
	E1 = C_UxMean[0]*X**2 + C_UxMean[1]*X + C_UxMean[2]
	E2 = C_UyMean[0]*X**2 + C_UyMean[1]*X + C_UyMean[2]
	E3 = C_uxRMS[0]*X**2 + C_uxRMS[1]*X + C_uxRMS[2]
	E4 = C_uyRMS[0]*X**2 + C_uyRMS[1]*X + C_uyRMS[2]
	E5 = C_uv[0]*X**2 + C_uv[1]*X + C_uv[2]
	mpl.plot(X,E1)
#	mpl.plot(z,error_UxMean)
	mpl.show()
	mpl.plot(X,E2)
#	mpl.plot(z,error_UyMean)
	mpl.show()
	mpl.plot(X,E3)
	mpl.plot(z,error_uxRMS)
	mpl.show()
	mpl.plot(X,E4)
	mpl.plot(z,error_uyRMS)
	mpl.show()
	mpl.plot(X,E5)
	mpl.plot(z,error_uv)
	mpl.show()
##
##
#########	New function: errorCompiler
##		Inputs are a list of data frames and a write name
##		Function outputs a new dataframe consisting of Z location and errors for each variable.
##		Write name will consist of pump speed and X location
def errorCompiler(df_names,write_name):
#	Need to loop through df_names and append errors and locations
#	initialise:
	z,e_UxMean,e_UyMean,e_uxRMS,e_uyRMS,e_uv = [],[],[],[],[],[]
	for i in range(len(df_names)):
		data = pd.read_pickle(df_names[i])
#		print(data)
		z.append(data.NXYZ[3])
		e_UxMean.append(data.error_UxMean[0])
		e_UyMean.append(data.error_UyMean[0])
		e_uxRMS.append(data.error_uxRMS[0])
		e_uyRMS.append(data.error_uyRMS[0])
		e_uv.append(data.error_uv[0])
#
##	Now compile into a new dataFrame
	z = pd.Series(z)
	e_UxMean = pd.Series(e_UxMean)
	e_UyMean = pd.Series(e_UyMean)
	e_uxRMS = pd.Series(e_uxRMS)
	e_uyRMS = pd.Series(e_uyRMS)
	e_uv = pd.Series(e_uv)
#
	errorData = pd.DataFrame({'z':z,'error_UxMean':e_UxMean,'error_UyMean':e_UyMean,'error_uxRMS':e_uxRMS,'error_uyRMS':e_uyRMS,'error_uv':e_uv})
#	print(errorData)
	errorData.to_pickle(write_name)
#	
#	Now we write this to the location provided
	return errorData

def simplePlotter(**kargs):
	u1 = kargs['u1'];	u2 = kargs['u2'];	u3 = kargs['u3'];	u4 = kargs['u4'];	u5 = kargs['u5'];	u6 = kargs['u6'];
	v1 = kargs['v1'];	v2 = kargs['v2'];	v3 = kargs['v3'];	v4 = kargs['v4'];	v5 = kargs['v5'];	v6 = kargs['v6'];
	u1Lab = kargs['u1Lab'];	u2Lab = kargs['u2Lab'];	u3Lab = kargs['u3Lab'];	u4Lab = kargs['u4Lab'];	u5Lab = kargs['u5Lab'];	u6Lab = kargs['u6Lab'];
	ylabel = kargs['ylabel'];	xlabel = kargs['xlabel'];
	legend = kargs['legend'];
	writeString=kargs['writeString'];
#
	mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	mpl.rc('text', usetex=True)
	fig = mpl.figure()
	ax1 = fig.add_subplot(1,1,1)
	plot1, = ax1.plot(v1,u1,color='k',linestyle='-',label=u1Lab)#,linewidth='2')
	if isinstance(u2,pd.Series): 
		plot2, = ax1.plot(v2,u2,color='b',linestyle='-',label=u2Lab)
	if isinstance(u3,pd.Series): 
		plot3, = ax1.plot(v3,u3,color='g',linestyle='-',label=u3Lab)
	if isinstance(u4,pd.Series): 
		plot4, = ax1.plot(v4,u4,color='r',linestyle='-',label=u4Lab)
	if isinstance(u5,pd.Series):
		plot5, = ax1.plot(v5,u5,color='m',linestyle='-',label=u5Lab)
	if isinstance(u6,pd.Series): 
		plot6, = ax1.plot(v6,u6,color='c',linestyle='-',label=u6Lab)
	ax1.set_xlabel(xlabel,fontsize='30')
	ax1.set_ylabel(ylabel,fontsize='35',rotation=0,labelpad=45)
#	h.set_rotation(0)
	if legend == True:
		mpl.legend(handles=[plot1, plot2, plot3, plot4, plot5, plot6],loc=4,prop={'size':18},ncol=1)
#
	mpl.axis([0, np.max(v6),0.7, 1.3])
	mpl.minorticks_on()
	mpl.grid(True, which='minor',alpha=0.6)
	mpl.grid(True, which='major',linewidth=0.9)
	writePath = writeString
	mpl.savefig(writePath)
	mpl.show()
	mpl.close()



###########################################################################################
