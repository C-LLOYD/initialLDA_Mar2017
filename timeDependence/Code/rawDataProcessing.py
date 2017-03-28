###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the time-dependence tests of Mar2017.
####		Data is currently processed using two functions:
####
####			txtToDataFrame:	This function reads in the txt files and returns a
####					data frame consisting of probe position, sample
####					number, time stamp and two velocity components.
####
####			rawToProcessed:	Adds additional columns onto the dataFrame such as
####					means, RMS velocities and reynolds stresses.
####
####		Currently working on the plotting scripts. The modularity of the script
####		will later allow the addition of filtering operations to the raw time 
####		series.
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
####			column 4:	Ux
####			column 5:	Uy
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
#			column 2:	timeInterval (ignore this)
#			column 3:	Ux
#			column 4:	Uy
#
	sampleNumber,timeStamp,Ux,Uy = [],[],[],[]
	for i in range(index,len(content)):
		sampleNumber.append(content[i].split('\t')[0])
		timeStamp.append(content[i].split('\t')[1])
		Ux.append(content[i].split('\t')[3])
		Uy.append(content[i].split('\t')[4])
#
#	Store variables as lists of floats
	sampleNumber=[float(i) for i in sampleNumber]
	timeStamp=[float(i) for i in timeStamp]
#	change units to seconds
	timeStamp = np.divide(timeStamp,1000)
	Ux=[float(i) for i in Ux]
	Uy=[float(i) for i in Uy]
#
#	Change data to a dataFrame:
#	First change each list to a series, then combine the series.
#	Output the final data frame
	NXYZ = pd.Series(NXYZ)
	sampleNumber = pd.Series(sampleNumber)
	timeStamp = pd.Series(timeStamp)
	Ux = pd.Series(Ux)
	Uy = pd.Series(Uy)
	data = pd.DataFrame({'NXYZ':NXYZ,'sampleNumber':sampleNumber,'timeStamp':timeStamp,'Ux':Ux,'Uy':Uy})
	return data;

###########################################################################################
####	Function Definition:		rawToProcessed				###########
####
####	This function calculates statistical quantities from the data frame input (as 
####	described above). The statistical quantities are added as additional columns to the
####	original data frame.
####
####	Statistical quantities are calculated for each data entry using the previous values
####	i.e this is used to check convergence of stats.
####
####		Current stats: UxMean, UyMean, uxRMS, uyRMS, uv
####
###########################################################################################
def rawToProcessed (data):
#
##	Read in data frame and convert required variables to lists: This speeds up looping
##	process.
	Ux = data.Ux.tolist()
	Uy = data.Uy.tolist()
#
##	Initialise statistics and loop through the length of the time series
	UxMean,UyMean,uxRMS,uyRMS,uv = [],[],[],[],[]
	for N in range(len(Ux)):
		print(N, 'of', len(Ux))
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
		uxRMS.append(np.divide(sum(np.power(Ux-UxMeanNew,2)),N+1))
		uyRMS.append(np.divide(sum(np.power(Uy-UyMeanNew,2)),N+1))
		uv.append(np.divide(sum((Ux-UxMeanNew)*(Uy-UyMeanNew)),N+1))
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
	return data;

###########################################################################################

fileName = "../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000001.txt"
data = txtToDataFrame(fileName)
data = rawToProcessed(data)


#	Add some lines at +/- 5%
#plus5=[np.mean(Ux)*0.05+np.mean(Ux)]*len(Ux)
#min5 = [np.mean(Ux)*-0.05+np.mean(Ux)]*len(Ux)
#plus1=[np.mean(Ux)*0.01+np.mean(Ux)]*len(Ux)
#min1 = [np.mean(Ux)*-0.01+np.mean(Ux)]*len(Ux)
#mpl.plot(timeStamp,(np.divide(UxMean-np.mean(Ux),np.mean(Ux))))
#mpl.plot(timeStamp,uxRMS)
#mpl.plot(timeStamp,plus5)
#mpl.plot(timeStamp,min5)
#mpl.plot(timeStamp,plus1)
#mpl.plot(timeStamp,min1)
#mpl.xlabel('t')
#mpl.ylabel('mean(Ux)')
#mpl.savefig('/usr/not-backed-up/convergence.png')
#mpl.show()

#######################################################################################
####	Now run the functions


