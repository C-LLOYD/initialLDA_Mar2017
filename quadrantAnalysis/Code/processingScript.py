###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the time-dependence tests of Mar2017.
####		Data is currently processed using two functions:
####
####			txtToDataFrame:	This function reads in the txt files and returns a
####					data frame consisting of probe position, sample
####					number, time stamp, residence time, and two
####					velocity components.
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
#import re
import matplotlib.pyplot as plt
import pandas as pd
from processingFunctions import txtToDataFrame
from processingFunctions import rawToProcessed_weighted
from FilterFunctions import Filter

def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return array[idx]


#data=pd.read_pickle("../Data/processedData/dataFrames/4hz_x_400_z_1_data_filtered_moving_average_weighted.pkl")
plotTimeDep = False
plotFilDep = False
writeDataFrames = False
readData = True

if writeDataFrames == True:
	fileName = "../Data/rawData/Run9_x400_fl_8_hz_300secs.000001.txt"
	data_raw = txtToDataFrame(fileName,"../Data/processedData/")
	data_fil_ma = Filter(data_raw,'movingAverageFilter','mean',10,False,"../Data/processedData/")
	data_fil_ma_wei 	= 	rawToProcessed_weighted(data_fil_ma,False,"../Data/processedData/",'_data_filtered_moving_average_weighted.pkl')

	print(data_fil_ma)

if readData == True:
#	data = pd.read_pickle("../Data/processedData/x_400_z_1_data_filtered_moving_average_weighted.pkl")
	data = pd.read_pickle("../../timeDependence/Data/processedData/dataFrames/8hz_x_400_z_15_data_filtered_moving_average_weighted.pkl")
#
##	For now, take u component of velocity
U = data.Ux.as_matrix()
V = data.Uy.as_matrix()
tau = data.resTime.as_matrix()
u = ((U-np.sum(U*tau)/np.sum(tau)))#*tau/np.max(tau)
v = ((V-np.sum(V*tau)/np.sum(tau)))#*tau/np.max(tau)
uv = u*v

pdfTests = False
if pdfTests == True:
	#N = data.sampleNumber.as_matrix()
	N = len(u)
	#
	##	Set up bins such that they cover whole range of u
	Nbins = 49
	Delta = np.linspace(np.min(u),np.max(u),Nbins+1)
	DeltaCentre = (Delta[0:-1]+Delta[1:])/2
	n = np.zeros(len(Delta)-1)
	taui = np.zeros(len(Delta)-1)
	#
	for i in range(len(n)-1):
		binMin = Delta[i]; binMax = Delta[i+1];
		index = map(lambda x: binMin <= x < binMax, u)
		n[i] = np.sum(list(index))
		if n[i] > 0:
			taui[i] = np.mean(tau[list(index)])
#
	binEnd = (Delta[-2],Delta[-1])
	index = map(lambda x: binEnd[0] <= x < binEnd[1], u)
	n[-1] = np.sum(list(index))
	if n[i] > 0:
		taui[-1] = np.mean(tau[list(index)])
#
	PDF = n*taui/((Delta[1]-Delta[0])*(np.sum(taui)))
	PDF = n/((Delta[1]-Delta[0])*N)
	PDF = n/(N)
	#plt.plot(DeltaCentre/np.max(tau),PDF)
	plt.plot(DeltaCentre,n)
	#plt.hist(u/np.max(tau),bins=99,density=True)
	plt.show()
	plt.close()
	#[pdf,edges] = np.histogram(U*tau/np.max(tau),bins=99,normed=True)
	plt.hist(U,bins=49)#,normed=True)
	plt.show()
	plt.close()
	
corrTests = False
if corrTests == True:
#	print(data)
	xMax = np.max(u)
	xMin = -xMax
	yMax = xMax
	yMin = -yMax
#
##	define grid
	xEdg = np.linspace(xMin,xMax,50)
	xCen = (xEdg[0:-1]+xEdg[1:])/2	
	yEdg = np.linspace(yMin,yMax,50)
	yCen = (yEdg[0:-1]+yEdg[1:])/2
#	pdf = np.zeros()	
	XEdg,YEdg = np.meshgrid(xEdg,yEdg)
	XCen,YCen = np.meshgrid(xCen,yCen)
	PDF = np.zeros((len(XCen),len(YCen)))
	test= np.zeros((len(XCen),len(YCen)))
#
##	Set up hole tests - for each hole level, calculate 	
#
##	
	plt.scatter(u,v)
#	plt.show()
	for i in range(len(u)):
		#
		##	Extract coordinate
		coord = [u[i],v[i]]
		t = tau[i]
		A = find_nearest(xCen,coord[0])
		B = find_nearest(yCen,coord[1])
		Aloc = np.where(xCen == A)
		Bloc = np.where(yCen == B)
		PDF[Bloc,Aloc] = PDF[Bloc,Aloc] + t
		test[Bloc,Aloc] = test[Bloc,Aloc] + 1
#		print(t,test[Bloc,Aloc],PDF[Bloc,Aloc])
	PDF = (np.ma.masked_where(PDF==0,PDF))/sum(tau)
	test = (np.ma.masked_where(test==0,test))/len(u)
	WPDF = (PDF*XCen*YCen)
	test2 = test*XCen*YCen
#
	plt.figure()
#	fig1 = plt.pcolormesh(XCen,YCen,PDF)
	fig1 = plt.contour(XCen,YCen,abs((PDF-test)/PDF),10)
	plt.clabel(fig1, inline=1, fontsize=10)
	plt.axvline();	plt.axhline();

	plt.figure()
#	fig1 = plt.pcolormesh(XCen,YCen,PDF)
	fig3 = plt.contour(XCen,YCen,test,10)
#	plt.clabel(fig1, inline=1, fontsize=10)
	plt.axvline();	plt.axhline();

	plt.figure()
#	fig1 = plt.pcolormesh(XCen,YCen,PDF)
	fig2 = plt.contour(XCen,YCen,WPDF,10)
#	plt.clabel(fig2, inline=1, fontsize=10)
	plt.axvline();	plt.axhline();

	plt.figure()
#	fig1 = plt.pcolormesh(XCen,YCen,PDF)
	fig2 = plt.contour(XCen,YCen,test2,10)
#	plt.clabel(fig2, inline=1, fontsize=10)
	plt.axvline();	plt.axhline();

	plt.figure()
	plt.scatter(XCen,PDF)
	plt.scatter(XCen,test)
	plt.show()
	

##	Set up hole tests - for each hole level, calculate the contributions from each
#					quadrant, and time spent inside the hole
holeTests = True
if holeTests == True:
#	print(data)
	xMax = np.max(u)
	xMin = -xMax
	yMax = xMax
	yMin = -yMax
#
##	define grid
	xEdg = np.linspace(xMin,xMax,50)
	xCen = (xEdg[0:-1]+xEdg[1:])/2	
	yEdg = np.linspace(yMin,yMax,50)
	yCen = (yEdg[0:-1]+yEdg[1:])/2
#	pdf = np.zeros()	
	XEdg,YEdg = np.meshgrid(xEdg,yEdg)
	XCen,YCen = np.meshgrid(xCen,yCen)
	PDF = np.zeros((len(XCen),len(YCen)))
	WPDF = np.zeros((len(XCen),len(YCen)))
#
##	initialise H, t and quadrants
	H = np.linspace(0,20,1)
	print(H)
	Q1 = np.zeros(len(H))
	Q2 = np.zeros(len(H))
	Q3 = np.zeros(len(H))
	Q4 = np.zeros(len(H))
	inner = np.zeros(len(H))
	time = np.zeros(len(H))
##	Loop through length of u/v vectors and sort them into bins
	for k in range(len(H)):
		for i in range(len(u)):
		#
		##	Extract coordinate
			coord = [u[i],v[i]]#
			t = tau[i]
			A = find_nearest(xCen,coord[0])
			B = find_nearest(yCen,coord[1])
			Aloc = np.where(xCen == A)
			Bloc = np.where(yCen == B)
			PDF[Bloc,Aloc] = PDF[Bloc,Aloc] + t
#			WPDF[Bloc,Aloc] = PDF[Bloc,Aloc] + t*u[i]*v[i]
		PDF = (np.ma.masked_where(PDF==0,PDF))/sum(tau)
#		WPDF = (np.ma.masked_where(WPDF==0,WPDF))/sum(tau)
		WPDF = (PDF*XCen*YCen)
		print(np.sum(WPDF),np.sum(u*v*tau)/np.sum(tau),data.uv.as_matrix()[-1])
		#
		##	FIND THE HOLE
		
#	


#	
#
#######################################################################################



