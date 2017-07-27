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
#print(data)

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
	PDF = (np.ma.masked_where(PDF==0,PDF))/sum(tau)
	WPDF = (PDF*XCen*YCen)
#
	plt.figure()
#	fig1 = plt.pcolormesh(XCen,YCen,PDF)
	fig2 = plt.contour(XCen,YCen,WPDF,10)
#	plt.clabel(fig2, inline=1, fontsize=10)
	plt.axvline();	plt.axhline();
#
##	Set up hole tests - for each hole level, calculate the contributions from each
#					quadrant, and time spent inside the hole
holeTests = False
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
	H = np.linspace(0,20,21)
	Q1 = np.zeros(len(H))
	Q2 = np.zeros(len(H))
	Q3 = np.zeros(len(H))
	Q4 = np.zeros(len(H))
	inner = np.zeros(len(H))
	outer = np.zeros(len(H))
	time = np.zeros(len(H))
##	Loop through length of u/v vectors and sort them into bins
	for i in range(len(u)):
	#
	##	Extract coordinate and evaluate PDF
		coord = [u[i],v[i]]
		t = tau[i]
		A = find_nearest(xCen,coord[0])
		B = find_nearest(yCen,coord[1])
		Aloc = np.where(xCen == A)
		Bloc = np.where(yCen == B)
		PDF[Bloc,Aloc] = PDF[Bloc,Aloc] + t
	#PDF = (np.ma.masked_where(PDF==0,PDF))/sum(tau)
	PDF = PDF/sum(tau)
	WPDF = (PDF*XCen*YCen)
#		print(np.sum(WPDF),np.sum(u*v*tau)/np.sum(tau),data.uv.as_matrix()[-1])
		#
	for k in range(len(H)):
		##	Now split up the outer contributions into each quadrant
		innerPDF = np.where(np.abs(XCen*YCen) < np.abs(H[k]*np.sum(WPDF)), WPDF, 0)
		outerPDF = np.where(np.abs(XCen*YCen) >= np.abs(H[k]*np.sum(WPDF)), WPDF, 0)
		Q1PDF = np.where((XCen > 0) & (YCen > 0) & (np.abs(XCen*YCen) >= np.abs(H[k]*np.sum(WPDF))), WPDF, 0)
		Q2PDF = np.where((XCen < 0) & (YCen > 0) & (np.abs(XCen*YCen) >= np.abs(H[k]*np.sum(WPDF))), WPDF, 0)
		Q3PDF = np.where((XCen < 0) & (YCen < 0) & (np.abs(XCen*YCen) >= np.abs(H[k]*np.sum(WPDF))), WPDF, 0)
		Q4PDF= np.where((XCen > 0) & (YCen < 0) & (np.abs(XCen*YCen) >= np.abs(H[k]*np.sum(WPDF))), WPDF, 0)
#
##
#		plt.figure()
#		test = plt.contour(XCen,YCen,Q1PDF)
#		test2 = plt.contour(XCen,YCen,Q2PDF)
#		test3 = plt.contour(XCen,YCen,Q3PDF)
#		test4 = plt.contour(XCen,YCen,Q4PDF)
#		plt.show()
#
##
		inner[k] = np.sum(innerPDF)/np.sum(WPDF)
		outer[k] = np.sum(outerPDF)/np.sum(WPDF)
		Q1[k] = np.sum(Q1PDF)/np.sum(WPDF)
		Q2[k] = np.sum(Q2PDF)/np.sum(WPDF)
		Q3[k] = np.sum(Q3PDF)/np.sum(WPDF)
		Q4[k] = np.sum(Q4PDF)/np.sum(WPDF)
#
##	
	plt.figure()
	plot1, = plt.plot(H,Q1,'r',label='Q1')
	plot2, = plt.plot(H,Q2,'m',label='Q2')
	plot3, = plt.plot(H,Q3,'c',label='Q3')
	plot4, = plt.plot(H,Q4,'g',label='Q4')
	plot5, = plt.plot(H,inner,'b',label='Inner')
	plot6, = plt.plot(H,outer,'k',label='Outer')
	plt.legend(handles=[plot1, plot2, plot3, plot4, plot5, plot6],loc=5,prop={'size':12},ncol=2)
	plt.axvline();	plt.axhline();
#
##
	tempPDF = WPDF[:]
#	tempPDF[tempPDF == 0] = np.nan
	plt.figure()
	fig2 = plt.contourf(XCen,YCen,tempPDF,10)
	plt.axvline();	plt.axhline();
	plt.show()		
#	

#
##	NEW METHOD FOR IDENTIFYING QUADRANTS - WE DON'T USE THE PDF HERE
newQuadMethod = True
if newQuadMethod == True:
#	Quad = np.zeros(len(U))
#	Quad = np.where( (u > 0) & (v > 0), 1,Quad)
#	Quad = np.where( (u < 0) & (v > 0), 2,Quad)
#	Quad = np.where( (u < 0) & (v < 0), 3,Quad)
#	Quad = np.where( (u > 0) & (v < 0), 4,Quad)
	H = np.linspace(0,20,21)
	Q1H = np.zeros(len(H))
	Q2H = np.zeros(len(H))
	Q3H = np.zeros(len(H))
	Q4H = np.zeros(len(H))
	innerH = np.zeros(len(H))
	outerH = np.zeros(len(H))
	Q1Hres = np.zeros(len(H))
	Q2Hres = np.zeros(len(H))
	Q3Hres = np.zeros(len(H))
	Q4Hres = np.zeros(len(H))
	innerHres = np.zeros(len(H))
	outerHres = np.zeros(len(H))
	for k in range(len(H)):
#
##		FIND INDICES FOR EACH VARIABLE	-	Currently this section
#		Does not work - the commented variables below ARE correct but not properly implemented into the
#		H vector structure. Next work will rearrange this code. 
		Q1i = (u > 0) & (v > 0) & (np.abs(u*v*tau/np.mean(tau)) >= np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
		Q2i = (u < 0) & (v > 0) & (np.abs(u*v*tau/np.mean(tau)) >= np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
		Q3i = (u < 0) & (v < 0) & (np.abs(u*v*tau/np.mean(tau)) >= np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
		Q4i = (u > 0) & (v < 0) & (np.abs(u*v*tau/np.mean(tau)) >= np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
		inneri = (np.abs(u*v*tau/np.mean(tau)) < np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
		outeri = (np.abs(u*v*tau/np.mean(tau)) >= np.abs(H[k]*np.sum(u*v*tau)/np.sum(tau)))
#
##		NOW UPDATE NEW VARIABLES
		Q1H[k] = (np.sum((u*v*tau)[Q1i])/np.sum(tau[Q1i]))/(np.sum(u*v*tau)/np.sum(tau))
		Q2H[k] = (np.sum((u*v*tau)[Q2i])/np.sum(tau[Q2i]))/(np.sum(u*v*tau)/np.sum(tau))
		Q3H[k] = (np.sum((u*v*tau)[Q3i])/np.sum(tau[Q3i]))/(np.sum(u*v*tau)/np.sum(tau))
		Q4H[k] = (np.sum((u*v*tau)[Q4i])/np.sum(tau[Q4i]))/(np.sum(u*v*tau)/np.sum(tau))
		innerH[k] = (np.sum((u*v*tau)[inneri])/np.sum(tau[inneri]))/(np.sum(u*v*tau)/np.sum(tau))
		outerH[k] = (np.sum((u*v*tau)[outeri])/np.sum(tau[outeri]))/(np.sum(u*v*tau)/np.sum(tau))
		Q1Hres[k] = np.sum(tau[Q1i])/np.sum(tau)
		Q2Hres[k] = np.sum(tau[Q2i])/np.sum(tau)
		Q3Hres[k] = np.sum(tau[Q3i])/np.sum(tau)
		Q4Hres[k] = np.sum(tau[Q4i])/np.sum(tau)
		innerHres[k] = np.sum(tau[inneri])/np.sum(tau)
		outerHres[k] = np.sum(tau[outeri])/np.sum(tau)
	print(Q1H[2]+Q2H[2]+Q3H[2]+Q4H[2],outerH[2])
	print((Q1H[2]+Q2H[2]+Q3H[2]+Q4H[2]+innerH[2]))
#	print(Q1Hres+Q2Hres+Q3Hres+Q4Hres+innerHres)
#	print(outerHres+innerHres)
	plt.figure()
	plot1, = plt.plot(H,Q1H,'r',label='Q1')
	plot2, = plt.plot(H,Q2H,'m',label='Q2')
	plot3, = plt.plot(H,Q3H,'c',label='Q3')
	plot4, = plt.plot(H,Q4H,'g',label='Q4')
	plot5, = plt.plot(H,innerH,'b',label='Inner')
	plot6, = plt.plot(H,outerH,'k',label='Outer')
	plt.legend(handles=[plot1, plot2, plot3, plot4, plot5, plot6],loc=5,prop={'size':12},ncol=2)
	plt.axvline();	plt.axhline();
	plt.show()



#		Q1 = np.empty(len(u))*np.nan
#		Q2, Q3, Q4 = np.copy(Q1), np.copy(Q1), np.copy(Q1)
#		Q1res, Q2res, Q3res, Q4res =  np.copy(Q1), np.copy(Q1), np.copy(Q1), np.copy(Q1)
#		Q1[Q1i] = (u*v*tau)[Q1i]
#	print(Q1)
#		Q2[Q2i] = (u*v*tau)[Q2i]
##	print(Q2)
#		Q3[Q3i] = (u*v*tau)[Q3i]
#		Q4[Q4i] = (u*v*tau)[Q4i]
#		Q1res[Q1i] = tau[Q1i]
#		Q2res[Q2i] = tau[Q2i]
#		Q3res[Q3i] = tau[Q3i]
#		Q4res[Q4i] = tau[Q4i]
#	print(np.nansum(Q1res)/np.sum(tau),np.nansum(Q2res)/np.sum(tau),np.nansum(Q3res)/np.sum(tau),np.nansum(Q4res)/np.sum(tau))
#	print((np.nansum(Q1)+np.nansum(Q2)+np.nansum(Q3)+np.nansum(Q4))/np.sum(tau),np.sum(u*v*tau)/np.sum(tau))
#
##	Quadrants have been identified - do we worry about the hole now ...?

#	
#
#######################################################################################



