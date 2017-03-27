###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the time-dependence tests of Mar2017.
####		Data is read in, processed, and saved as separate data files for plotting.
####
###########################################################################################
####
##		Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
#
##		Read in data (This will later be looped as a function or otherwise)
#fileName = "../Data/rawData/4hz/Run8_x400_fl_4hz_300secs/Run8_x400_fl_4_hz_300secs.000003.txt"
def txtToArray (fileName):
	with open(fileName) as f:
		content=f.readlines()
		content = [x.strip() for x in content]
##		Content is a list of strings, separated by \n, containing each lines
##		First few lines represent header
##		Need to identify XYZ from header
##		Need to identify row number containing first data point
##
##		First grab the file number and position from the content list
	i=0
	while len(content[i].split(';')) is not 4:	#Test to see if the current line, when split by ';', has a length of 4
		i=i+1					#A length of 4 will correspond to the line containing NXYZ
	NXYZ = content[i].split(';')			#Stores the row containing the position and splits it by ';'
	for i in range(len(NXYZ)):			#loop through the list NXYZ	
		NXYZ[i] = re.sub(' mm$','',NXYZ[i])	#remove units and store only scalars
#		NXYZ is a list containing only numbers, corresponding to file number and XYZ positions.
##
##		Now check to find the raw data
##		Data column names are split by tabs and the first entry is always "Row#"
##		We use another while loop to find this row number and then we add one to find the data
	i=0
	for i in range(len(content)):
		if re.match(r'"Row#"',content[i].split('\t')[0]):
			index=i+1
			break
##		Now we have located this start of the data by the index variable
##		Now loop through the length of content and extract the important variables
##
##		DATA ARRANGED AS:	sampleNumber,	timeStamp,	timeInterval,	Ux,	Uy
##		We don't need to extract timeInterval but the others are useful.
##
#print(content[index])
	sampleNumber,timeStamp,Ux,Uy = [],[],[],[]
	for i in range(index,len(content)):
		sampleNumber.append(content[i].split('\t')[0])
		timeStamp.append(content[i].split('\t')[1])
		Ux.append(content[i].split('\t')[3])
		Uy.append(content[i].split('\t')[4])
	sampleNumber=[float(i) for i in sampleNumber]
	timeStamp=[float(i) for i in timeStamp]
	timeStamp = np.divide(timeStamp,1000)		#change units to seconds
	Ux=[float(i) for i in Ux]
	Uy=[float(i) for i in Uy]
	dataStack = np.vstack([sampleNumber,timeStamp,Ux,Uy])
	return	dataStack;


########################################################################################
####
####	We now have 4 vectors containing raw data - now we create new variables
##
##	To compute means and stresses etc. we loop through the number of points and
##	compute quantities based on the previous entries i.e check convergence as we 
##	include more points / timeStamps.

fileName = "../Data/rawData/8hz/Run9_x400_fl_8hz_300secs/Run9_x400_fl_8_hz_300secs.000001.txt"
dataStack = txtToArray(fileName)
sampleNumber=np.ndarray.tolist(dataStack[0])
timeStamp=np.ndarray.tolist(dataStack[1])
Ux=np.ndarray.tolist(dataStack[2])
Uy=np.ndarray.tolist(dataStack[3])

UxMean,UyMean,uxRMS,uyRMS,uv = [],[],[],[],[]
for N in range(0,len(sampleNumber)):
#	Compute means
#	
#	print N, 'of', len(sampleNumber)
	UxMeanNew = []
	UyMeanNew = []
	UxMeanNew = np.divide(sum(Ux[0:N+1]),(N+1))
	UxMean.append(UxMeanNew)
	UyMeanNew = np.divide(sum(Uy[0:N+1]),(N+1))
	UyMean.append(UyMeanNew)
#
	uxRMS.append(np.divide(sum(np.power(Ux-UxMeanNew,2)),N+1))

print np.sqrt(uxRMS[len(uxRMS)-1])
print np.std(Ux)

#	Add some lines at +/- 5%
plus5=[np.mean(Ux)*0.05+np.mean(Ux)]*len(Ux)
min5 = [np.mean(Ux)*-0.05+np.mean(Ux)]*len(Ux)
plus1=[np.mean(Ux)*0.01+np.mean(Ux)]*len(Ux)
min1 = [np.mean(Ux)*-0.01+np.mean(Ux)]*len(Ux)
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


