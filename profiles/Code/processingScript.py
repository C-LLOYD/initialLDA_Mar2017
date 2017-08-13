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
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import glob
from processingFunctions import txtToDataFrame
from processingFunctions import timeAverage
from processingFunctions import plotter
from processingFunctions import findDimensionlessParameters as FDP
from FilterFunctions import Filter
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
writeData = False
#writeData = True
plotDimensionedData = False


rawPath = 	["../Data/rawData/4hz/300mm/*/*.txt",
			"../Data/rawData/4hz/400mm/*/*.txt",
			"../Data/rawData/8hz/300mm/*/*.txt",
			"../Data/rawData/8hz/400mm/*/*.txt",
			"../Data/rawData/16hz/300mm/*/*.txt",
			"../Data/rawData/16hz/400mm/*/*.txt"]
#
errorPath = 	["../Data/processedData/dataFrames/4hz_errors.pkl",
			"../Data/processedData/dataFrames/4hz_errors.pkl",
			"../Data/processedData/dataFrames/8hz_errors.pkl",
			"../Data/processedData/dataFrames/8hz_errors.pkl",
			"noData",
			"noData"]
#
writePath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/16hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/16hz_400mm_profiles.pkl"]
#
##
dataPath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/16hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/16hz_400mm_profiles.pkl"]
#
#
#
if writeData == True:
	for j in range(len(dataPath)):
		data = []
		print(rawPath[j])
		for fileName in glob.glob(rawPath[j]):
			tempData = txtToDataFrame(fileName)
#
##	3.	filter data
			if isinstance(tempData,pd.DataFrame):
				tempData = Filter(tempData,'movingAverageFilter','mean',10,'none','none')
#
##	4.	apply averaging and append a final data series
				dataNew  = timeAverage(tempData)
				if not isinstance(data,pd.DataFrame):
					data = dataNew
				else:
					data=data.append(dataNew)
		data["error_UxMean"] = np.nan
		data["error_UyMean"] = np.nan
		data["error_uxRMS"] = np.nan
		data["error_uyRMS"] = np.nan
		data["error_uv"] = np.nan
		data = data.sort_values(by='z',ascending=1)
		data = data.set_index("z")
#
##	This section of code reads in the real error data and interpolates over missing entries. First we read in the error data, then
##	rearrange for the index 'z', then place the corresponding values into each error vector. Last step - interpolate over missing data.
		if errorPath[j] is not "noData": 
			errorData = pd.read_pickle(errorPath[j])
			errorData = errorData.set_index("z")
			exactZ = errorData.index
			for i in range(len(exactZ)):
				data.set_value(exactZ[i],"error_UxMean", 	errorData.error_UxMean.loc[errorData.index == exactZ[i]])
				data.set_value(exactZ[i],"error_UyMean", 	errorData.error_UyMean.loc[errorData.index == exactZ[i]])
				data.set_value(exactZ[i],"error_uxRMS",	 	errorData.error_uxRMS.loc[errorData.index == exactZ[i]])
				data.set_value(exactZ[i],"error_uyRMS", 	errorData.error_uyRMS.loc[errorData.index == exactZ[i]])
				data.set_value(exactZ[i],"error_uv", 		errorData.error_uv.loc[errorData.index == exactZ[i]])
#	
			data = data.reset_index(level=0)
			errorData = errorData.reset_index(level=0)
			data.set_value(0,"error_UxMean", 	errorData.error_UxMean.loc[errorData.z == np.min(errorData.z)])
			data.set_value(0,"error_UyMean", 	errorData.error_UyMean.loc[errorData.z == np.min(errorData.z)])
			data.set_value(0,"error_uxRMS",	 	errorData.error_uxRMS.loc[errorData.z == np.min(errorData.z)])
			data.set_value(0,"error_uyRMS", 	errorData.error_uyRMS.loc[errorData.z == np.min(errorData.z)])
			data.set_value(0,"error_uv", 		errorData.error_uv.loc[errorData.z == np.min(errorData.z)])
			data = data.set_index("z")
			data = data.apply(pd.Series.interpolate)
#
		data = data.reset_index(level=0)
#
		print(data)
		data.to_pickle(writePath[j])
#
#
##################################################################################################################################
####
####		Set up variables for plotting and fitting


plotHighSpeed = True
if plotHighSpeed == True:
	d1 = pd.read_pickle(dataPath[4])
	d2 = pd.read_pickle(dataPath[5])
##	Plot mean U
	f=plt.figure()
	ax = f.add_subplot(111)
	font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
	plt.rc('font', **font)
	plt.rc('text', usetex=True)
	plot1, = plt.semilogx(d1.z,d1.UxMean,linestyle = 'None',color = 'k', marker = '.',label = r'$x=300$ (mm)',markersize=8)
	plot2, = plt.semilogx(d2.z,d2.UxMean,linestyle = 'None',color = 'r', marker = 'x',label = r'$x=400$ (mm)',markersize=8)
	plt.xlabel(r'$y$ (mm)',fontsize=30)
	plt.ylabel(r'$\mu_u$ (m/s)',fontsize=30)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
	ax.tick_params(axis='x', pad=15)
	#	plt.legend(handles=[plot1,plot2],loc=2)
	write = '../Data/processedData/figures/meanU_16hz.png'
#		plt.savefig(write)
	plt.show()
	plt.close()
##	Plot std U
	f=plt.figure()
	ax = f.add_subplot(111)
	font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
	plt.rc('font', **font)
	plt.rc('text', usetex=True)
	plot1, = plt.semilogx(d1.z,d1.uxRMS,linestyle = 'None',color = 'r', marker = '.',label = r'$\sigma_u$ at $x=300$ (mm)',markersize=8)
	plot2, = plt.semilogx(d1.z,d1.uyRMS,linestyle = 'None',color = 'g', marker = '.',label = r'$\sigma_v$ at $x=300$ (mm)',markersize=8)
	#
	plot3, = plt.semilogx(d2.z,d2.uxRMS,linestyle = 'None',color = 'b', marker = 'x',label = r'$\sigma_u$ at $x=400$ (mm)',markersize=8)
	plot4, = plt.semilogx(d2.z,d2.uyRMS,linestyle = 'None',color = 'k', marker = 'x',label = r'$\sigma_v$ at $x=400$ (mm)',markersize=8)
	plt.xlabel(r'$y$ (mm)',fontsize=30)
	plt.ylabel(r'$\sigma_u,\sigma_v$ (m/s)',fontsize=30)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
	ax.tick_params(axis='x', pad=15)
	ax.set_ylim(bottom=0)
	#	plt.legend(handles=[plot1,plot2,plot3,plot4],loc=1)
	write = '../Data/processedData/figures/std_16hz.png'
#		plt.savefig(write)
	plt.show()
	plt.close()
##	Plot cov
	f=plt.figure()
	ax = f.add_subplot(111)
	font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
	plt.rc('font', **font)
	plt.rc('text', usetex=True)
	plot1, = plt.semilogx(d1.z,d1.uv,linestyle = 'None',color = 'k', marker = '.',label = r'$x=300$ (mm)',markersize=8)
	#
	plot2, = plt.semilogx(d2.z,d2.uv,linestyle = 'None',color = 'r', marker = 'x',label = r'$x=400$ (mm)',markersize=8)
	plt.xlabel(r'$y$ (mm)',fontsize=30)
	plt.ylabel(r'$\gamma_{u,v}$ (m\textsuperscript{2}/s\textsuperscript{2})',fontsize=30)
	ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
	ax.tick_params(axis='x', pad=15)
	mpl.ticker.ScalarFormatter(useMathText = True)
	write = '../Data/processedData/figures/cov_16hz.png'
	plt.savefig(write)
	plt.show()
#		plt.close()

##
	
#
##	Current code section only plots data for lower pump speeds - 16Hz speed is done separately since there is no error data
if plotDimensionedData == True:
	data = (pd.read_pickle(dataPath[0]),pd.read_pickle(dataPath[1]),pd.read_pickle(dataPath[2]),pd.read_pickle(dataPath[3]))
	data1=data[0];	data2=data[1];	data3=data[2];	data4=data[3];
	pumpSpeed = ['4hz','8hz']
	for i in range(len(pumpSpeed)):
		if i == 0:
			d1 = data1;	d2 = data2;
		elif i == 1:
			d1 = data3; d2 = data4;
#
##	Plot mean U
		f=plt.figure()
		ax = f.add_subplot(111)
		font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
		plt.rc('font', **font)
		plt.rc('text', usetex=True)
	#	plt.rc('text', usetex=True)
	#	plt.rc('font', family='serif')
		plot1, = plt.semilogx(d1.z,d1.UxMean,linestyle = 'None',color = 'k', marker = '.',label = r'$x=300$ (mm)',markersize=8)
		plt.errorbar(d1.z,d1.UxMean,yerr=d1.error_UxMean, linestyle = "None",color = 'k')
	#
		plot2, = plt.semilogx(d2.z,d2.UxMean,linestyle = 'None',color = 'r', marker = 'x',label = r'$x=400$ (mm)',markersize=8)
		plt.errorbar(d2.z,d2.UxMean,yerr=d2.error_UxMean, linestyle = "None",color = 'k')
	#	plt.xticks(fontsize=25)
	#	plt.yticks(fontsize=25)
		plt.xlabel(r'$y$ (mm)',fontsize=30)
		plt.ylabel(r'$\mu_u$ (m/s)',fontsize=30)
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
		ax.tick_params(axis='x', pad=15)
	#	plt.legend(handles=[plot1,plot2],loc=2)
		write = '../Data/processedData/figures/meanU_'+str(pumpSpeed[i])+'.png'
		plt.savefig(write)
#		plt.show()
		plt.close()
##	Plot std U
		f=plt.figure()
		ax = f.add_subplot(111)
		font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
		plt.rc('font', **font)
		plt.rc('text', usetex=True)
#		plt.rc('text', usetex=True)
#		plt.rc('font', family='serif')
		plot1, = plt.semilogx(d1.z,d1.uxRMS,linestyle = 'None',color = 'r', marker = '.',label = r'$\sigma_u$ at $x=300$ (mm)',markersize=8)
		plt.errorbar(d1.z,d1.uxRMS,yerr=d1.error_uxRMS, linestyle = "None",color = 'k')
		plot2, = plt.semilogx(d1.z,d1.uyRMS,linestyle = 'None',color = 'g', marker = '.',label = r'$\sigma_v$ at $x=300$ (mm)',markersize=8)
		plt.errorbar(d1.z,d1.uyRMS,yerr=d1.error_uyRMS, linestyle = "None",color = 'k')
	#
		plot3, = plt.semilogx(d2.z,d2.uxRMS,linestyle = 'None',color = 'b', marker = 'x',label = r'$\sigma_u$ at $x=400$ (mm)',markersize=8)
		plt.errorbar(d2.z,d2.uxRMS,yerr=d2.error_uxRMS, linestyle = "None",color = 'k')
		plot4, = plt.semilogx(d2.z,d2.uyRMS,linestyle = 'None',color = 'k', marker = 'x',label = r'$\sigma_v$ at $x=400$ (mm)',markersize=8)
		plt.errorbar(d2.z,d2.uyRMS,yerr=d2.error_uyRMS, linestyle = "None",color = 'k')
#
	#	plt.xticks(fontsize=25)
	#	plt.yticks(fontsize=25)
		plt.xlabel(r'$y$ (mm)',fontsize=30)
		plt.ylabel(r'$\sigma_u,\sigma_v$ (m/s)',fontsize=30)
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
		ax.tick_params(axis='x', pad=15)
		ax.set_ylim(bottom=0)
	#	plt.legend(handles=[plot1,plot2,plot3,plot4],loc=1)
		write = '../Data/processedData/figures/std_'+str(pumpSpeed[i])+'.png'
		plt.savefig(write)
#		plt.show()
		plt.close()
##	Plot cov
		f=plt.figure()
		ax = f.add_subplot(111)
		font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
		plt.rc('font', **font)
		plt.rc('text', usetex=True)
	#	plt.rc('font', family='serif')
	#	plt.rc('font', size=24)
	#	plt.rc('font', weight='normal')
		plot1, = plt.semilogx(d1.z,d1.uv,linestyle = 'None',color = 'k', marker = '.',label = r'$x=300$ (mm)',markersize=8)
		plt.errorbar(d1.z,d1.uv,yerr=d1.error_uv, linestyle = "None",color = 'k')
	#
		plot2, = plt.semilogx(d2.z,d2.uv,linestyle = 'None',color = 'r', marker = 'x',label = r'$x=400$ (mm)',markersize=8)
		plt.errorbar(d2.z,d2.uv,yerr=d2.error_uv, linestyle = "None",color = 'k')
	#	plt.xticks(fontsize=25)
	#	plt.yticks(fontsize=25)
		plt.xlabel(r'$y$ (mm)',fontsize=30)
		plt.ylabel(r'$\gamma_{u,v}$ (m\textsuperscript{2}/s\textsuperscript{2})',fontsize=30)
#		plt.legend(handles=[plot1,plot2],loc=3)
	#	ax.legend(handles=[plot1,plot2],loc='upper right', bbox_to_anchor=(1, 1.05))
	#	ax.tick_params(labelsize=24)
	#	plt.axis.get_major_formatter().set_powerlimits((0, 1))
		ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
		ax.tick_params(axis='x', pad=15)
		mpl.ticker.ScalarFormatter(useMathText = True)
		write = '../Data/processedData/figures/cov_'+str(pumpSpeed[i])+'.png'
		plt.savefig(write)
#		plt.show()
		plt.close()

################	TEMP LOCATION OF LAW OF THE WALL FUNCTION - MOVE THIS LATER ################
####
####	1. 	define new variables
##	Initialise variables using the full data set and iterate to find log-law region by removing extremes of z


def linearReg(X,Y):
	if len(X) < 5:
#		print('N < 5')
		return [0, 0, 0]
	else:
#
##	Equation: X = beta*Y + alpha
		beta = np.sum( (X-np.mean(X)) * (Y-np.mean(Y)) ) / np.sum( (Y - np.mean(Y))**2 ) 
		alpha = np.mean(X) - beta*np.mean(Y)
		r2 = ((np.mean(X*Y)-np.mean(X)*np.mean(Y))/np.sqrt((np.mean(X**2)-np.mean(X)**2)*(np.mean(Y**2)-np.mean(Y)**2)))**2
		return [alpha,beta,r2]

def wallFinder(U,Y):
##	Estimate linear region
		N = np.zeros(len(U))
		r = np.zeros(len(U))
		a = np.zeros(len(U))
		b = np.zeros(len(U))
#
##	Calculate r2 for whole range of N and locate maximum r2 value
##		This corresponds to the best fit
		for i in range(len(U)):
			[a0,b0,r0] = linearReg(Y[0:i],U[0:i])
			r[i] = r0
			a[i] = a0
			b[i] = b0
			N[i] = i
#		print(r,N)
#		plt.plot(N,r)
#		plt.show()
#		plt.close()
		amax = a[np.where(r == np.max(r))][0]
		bmax = b[np.where(r == np.max(r))][0]
		rmax = r[np.where(r == np.max(r))][0]
		Nmax = int(N[np.where(r == np.max(r))][0])
#		print(type(Nmax))
#
##	We now have laminar layer : U[0:Nmax] = bmax*Y[0:Nmax] + amax
##	
##	Use this to find the wall and subtract this from the y values
#		Delta = amax/bmax
		Ynew = Y - amax
		[a1,b1,r1] = linearReg(Ynew[0:Nmax],U[0:Nmax])
#
##	Estimate Utau from linear part
		nu = 1.1*10**(-6)
		Utau = np.sqrt(nu/b1)
#		plt.plot(Ynew[0:Nmax+10],U[0:Nmax+10],linestyle=' ',marker='o')
#		plt.plot(U[0:Nmax+10]*b1+a1,U[0:Nmax+10])
#		plt.show()
#		plt.close()
#		print(amax,a1,rmax,r1)
		return [Utau,Ynew]



test2 = False
if test2 == True:
	d = [pd.read_pickle(dataPath[0]),pd.read_pickle(dataPath[1]),pd.read_pickle(dataPath[2]),pd.read_pickle(dataPath[3])]
	case = ['4Hz_300','4Hz_400','8Hz_300','8Hz_400']
	for j in range(len(d)):
		U = d[j].UxMean.as_matrix()
		e = d[j].error_UxMean.as_matrix()
		Y = d[j].z.as_matrix()/1000
		[UtauEst,Ynew] = wallFinder(U,Y)
		Y = Ynew
#
##	Check if this works as a scalar:
		nu = 1.1*10**(-6)
		kappa = 0.41
		B = 5.2		
#		plt.semilogx(Y*UtauEst/nu,U/UtauEst,ls=' ',marker='o')
#		plt.plot(Y[0:100]*UtauEst/nu,Y[0:100]*UtauEst/nu)
#		plt.plot(Y*UtauEst/nu,np.log(Y*UtauEst/nu)/kappa+B)
#		plt.show()
#		plt.close()
#
##	Use as estimate for range of log-law
		Ulog = U[np.where((Y*UtauEst/nu > 30) & (Y*UtauEst/nu < 150))]
		Ylog = np.log(Y[np.where((Y*UtauEst/nu > 30) & (Y*UtauEst/nu < 150))])
#		print(Ulog,Ylog)
#		print(U,Y)
#		print(type(Ulog),type(U))
		[a,b,r] = linearReg(Ulog,Ylog)
		Utau = b*kappa
		print(Utau,UtauEst,(UtauEst-Utau)/UtauEst)#,0.1*Utau/nu,0.1*UtauEst/nu)
		B = a/Utau - np.log(Utau/nu)/kappa
#		print(B)
#
##		Plot new

		f=plt.figure()
		ax = f.add_subplot(111)
		font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
#
		plt.rc('font', **font)
		plt.rc('text', usetex=True)
#
		plt.semilogx(Y*Utau/nu,U/Utau,ls=' ',marker='o',color='k')
		plt.errorbar(Y*Utau/nu,U/Utau,yerr=e/Utau, linestyle = "None",color = 'k')
		plt.plot(Y[0:100]*UtauEst/nu,Y[0:100]*UtauEst/nu,ls='-',color='r')
		plt.plot(Y*Utau/nu,np.log(Y*Utau/nu)/kappa+B,ls='-',color='b')
		plt.xlabel(r'$y^+$',fontsize=30)
		plt.ylabel(r'$U^+$',fontsize=30)
#		ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,1))
		ax.tick_params(axis='x', pad=15)
		mpl.ticker.ScalarFormatter(useMathText = True)
		write = '../Data/processedData/figures/Uplus_Yplus_'+str(case[j])+'.png'
		plt.savefig(write)
#		plt.show()
		plt.close()
		

#		[Delta,Ynew] = wallFinder(U,Y)
#		print(Delta)
#		print(Ynew)
#


#
##
#	[a,b,r0] = linearReg(U[0:N],Y[0:N])
#	print(r0)
#	r1 = r0
#	while r1 <= r0:
#		N = N - 2
#		Ulin = U[0:N]
#		Ylin = Y[0:N]
#		r0 = r1
#		plt.plot(U[0:N],Y[0:N],ls=' ',marker='o')
#		plt.plot(Y[0:N]*b+a,Y[0:N],ls='-')	
#		[a,b,r1] = linearReg(U[0:N],Y[0:N])
#		plt.plot(Y[0:N]*b+a,Y[0:N],ls='--')
#		plt.show()
#		print(r1)
		
#
##		Shape factor and integral methods etc...
test3 = False
if test3 == True:
	d = [pd.read_pickle(dataPath[0]),pd.read_pickle(dataPath[1]),pd.read_pickle(dataPath[2]),pd.read_pickle(dataPath[3])]
	case = ['4Hz_300','4Hz_400','8Hz_300','8Hz_400']
	for j in range(len(d)):
		U = d[j].UxMean.as_matrix()
		Y = d[j].z.as_matrix()/1000
		U0 = np.max(U)
		deltaStar = np.sum((1-U/U0)*Y)
		theta = np.sum((1-U/U0)*U*Y/U0)
		H = deltaStar/theta
		ReTheta = theta*U0/1.1e-6
		print(case[j],ReTheta,H)
	




testing = False
if testing == True:
	data4 = pd.read_pickle(dataPath[0])
	data8 = pd.read_pickle(dataPath[2])

	Ulab4 = data4.UxMean.as_matrix()
	Ylab4 = data4.z.as_matrix()/1000
	Ulab8 = data8.UxMean.as_matrix()
	Ylab8 = data8.z.as_matrix()/1000
	
	nu = 1.1*10**(-6)
	kappa = 0.41
	Cplus = 5.2
	[Utau4,Delta4] = FDP(Ulab4,Ylab4)
	[Utau8,Delta8] = FDP(Ulab8,Ylab8)
	
	Yplus4 = (Ylab4-Delta4)*Utau4/nu
	Uplus4 = Ulab4/Utau4
	Yplus8 = (Ylab8-Delta8)*Utau8/nu
	Uplus8 = Ulab8/Utau8
	
	plt.semilogx(Yplus4,Uplus4,linestyle=' ',marker = 'o')
	plt.errorbar(Yplus4,Uplus4,yerr=data4.error_UxMean/Utau4, linestyle = "None",ecolor='k', capthick=2)
	#
	plt.plot(Yplus4,Yplus4)
	plt.plot(Yplus4,np.log(Yplus4)/kappa+Cplus)
	plt.axis([0,np.max(Yplus4),0,np.max(Uplus4)+0.05*np.max(Uplus4)])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)
	plt.xlabel(r'$\mathbf{z^+}$',fontsize=30)
	plt.ylabel(r'$\mathbf{U^+}$',fontsize=30)
	plt.savefig('../Data/processedData/figures/Uplus_4hz.png')
	#plt.show()
	plt.close()
	
	plt.semilogx(Yplus8,Uplus8,linestyle=' ',marker = '*')
	plt.errorbar(Yplus8,Uplus8,yerr=data8.error_UxMean/Utau8, linestyle = "None",ecolor='k', capthick=2)
	#
	plt.plot(Yplus8,Yplus8)
	plt.plot(Yplus8,np.log(Yplus8)/kappa+Cplus)
	plt.axis([0,np.max(Yplus8),0,np.max(Uplus8)+0.05*np.max(Uplus8)])
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	plt.xticks(fontsize=25)
	plt.yticks(fontsize=25)
	plt.xlabel(r'$\mathbf{z^+}$',fontsize=30)
	plt.ylabel(r'$\mathbf{U^+}$',fontsize=30)
	#plt.savefig('../Data/processedData/figures/Uplus_8hz.png')
	#plt.show()
	plt.close()
	
	
#	Test y values for a given y^+
#print('4hz: For yplus = 5: y = ',5*nu/Utau4)
#print('4hz: For yplus = 10: y = ',10*nu/Utau4)
#print('4hz: For yplus = 20: y = ',20*nu/Utau4)
#print('8hz: For yplus = 5: y = ',5*nu/Utau8)
#print('8hz: For yplus = 10: y = ',10*nu/Utau8)
#print('8hz: For yplus = 20: y = ',20*nu/Utau8)




checkUtau = False
if checkUtau == True:
	for i in np.linspace(Utau/2,Utau + Utau/2, num=20):
		Yplus = (Ylab-Delta)*i/nu
		Uplus = Ulab/i
		plt.semilogx(Yplus,Uplus,linestyle=' ',marker = '.')
		plt.plot(Yplus,Yplus)
		plt.plot(Yplus,np.log(Yplus)/kappa+Cplus)
		plt.axis([0,np.max(Yplus),0,np.max(Uplus)+0.05*np.max(Uplus)])
		plt.show()
		plt.close()


checkKappa = False
if checkKappa == True:
	for kappait in np.linspace(0.38,0.42,num=5):
		print(kappait)
		plt.semilogx(Yplus,Uplus,linestyle=' ',marker = '.')
		plt.plot(Yplus,Yplus)
		plt.plot(Yplus,np.log(Yplus)/kappait+Cplus)
		plt.axis([0,np.max(Yplus),0,np.max(Uplus)+0.05*np.max(Uplus)])
		plt.show()
		plt.close()

checkCplus = False
if checkCplus == True:
	for Cplusit in np.linspace(5,6,num=10):
		print(Cplusit)
		plt.semilogx(Yplus,Uplus,linestyle=' ',marker = '.')
		plt.plot(Yplus,Yplus)
		plt.plot(Yplus,np.log(Yplus)/kappa+Cplusit)
		plt.axis([0,np.max(Yplus),0,np.max(Uplus)+0.05*np.max(Uplus)])
		plt.show()
		plt.close()


#
##
#Ylam = Ylab[0:20]
#Ulam = Ylab[0:20]
#counter = 0
#while counter < 5:
#	counter = counter + 1
#	print(counter)
#	print(Ylam,Ulam)
#	[aLam,bLam] = linearReg(Ulam,Ylam)
#	Delta = aLam/bLam
#	Utau1 = np.sqrt(bLam*nu)
#	Delta = aLam*nu/(Utau1**2)
#	Yplus = (Ylab+Delta)*Utau1/nu
##
###	Use these to estimate range for log region
#	Yturb = Ylab[(30<=Yplus) & (Yplus<=200)]
#	Uturb = Ulab[(30<=Yplus) & (Yplus<=200)]
#
##	Now estimate Utau based on Delta and linear regression of log law component
#	[aturb,bTurb] = linearReg(Uturb,np.log(Yturb+Delta))
#	Utau2 = bTurb*kappa
#	Yplus = Ylab*Utau2/nu
#	Ylam = Ylab[(Yplus<=5)]
##	print(len(Ylam))
#	Ulam = Ulab[(Yplus<=5)]
#	if len(Ylam) < 2:
#		Ylam = Ylam[0:3]
#		Ulam = Ulam[0:3]
#	print(Utau1,Utau2)
#Utau = Utau2
#Yplus = (Ylab + Delta)*Utau/nu
#Uplus = Ulab/Utau
#
#lamlawUplus = Yplus
#loglawUplus = np.log(Yplus)/kappa + Cplus# + (Utau1/kappa)*np.log(Utau1/nu)#
#
#plt.semilogx(Yplus,Uplus,linestyle=' ',marker = '.')
#plt.scatter((Ylam+Delta)*Utau/nu,Ulam/Utau)
#plt.scatter((Yturb+Delta)*Utau/nu,Uturb/Utau)
#plt.semilogx(Yplus,lamlawUplus)
#plt.semilogx(Yplus,loglawUplus)
#plt.axis([0,np.max(Yplus),0,np.max(Uplus)+0.05*np.max(Uplus)])
#plt.show()
#plt.close()

#plt.plot(Uplus,Yplus,linestyle=' ',marker= '.')
#plt.plot(lamlawUplus,Yplus)
#plt.plot(loglawUplus,Yplus)	
#plt.axis([0,30,0,100])
#plt.show()
#plt.close()	

#
##
#U_lam = U[0:10]
#Y_lam = Y[0:10]
#
###	Compute alpha and beta for the laminar layer
#betaHat_lam = np.sum((U_lam-np.mean(U_lam))*(Y_lam-np.mean(Y_lam)))/np.sum((Y_lam-np.mean(Y_lam))**2)
#alphaHat_lam = np.mean(U_lam) - betaHat_lam*np.mean(Y_lam)
#
##	Use these to estimate offset from the wall and Utau
#nu=1.1*10**(-6)
#Delta = alphaHat_lam/betaHat_lam
#Utau = np.sqrt(betaHat_lam*nu)
#print(Utau)
#
##	Now remove the offset from Y
#Y = Y + Delta
#Yplus = Y*Utau/nu
#
##	Now we estimate which points lie in the log-law region
##		30 < yPlus_loglaw < 200
#Y_turb = Y[(30<=Yplus) & (Yplus<=200)]
#U_turb = U[(30<=Yplus) & (Yplus<=200)]
#betaHat_turb = np.sum((U_turb-np.mean(U_turb))*(np.log(Y_turb)-np.mean(np.log(Y_turb))))/np.sum((np.log(Y_turb)-np.mean(np.log(Y_turb)))**2)
#print(betaHat_turb)
#alphaHat_turb = np.mean(U_turb) - betaHat_turb*np.mean(np.log(Y_turb))
#
##
#ReTau_lam = Utau*0.1/nu
#ReTau_turb = (0.1/nu)*0.41*betaHat_turb
#
##	Assume free stream velocity is equal the the velocity at the last point
#Cf = 2*Utau**2/U[-1]**2
#
#print(ReTau_lam,ReTau_turb,Cf)
#
#Utau = betaHat_turb*0.41
#yPlus = (Y*Utau/nu)#/1000
#UPlus = U/Utau
#print(yPlus,UPlus)
#
##	Plot this
#loglawU = np.zeros(len(U))
#loglawU[:] = np.log(yPlus)/0.41 + 5
#			#betaHat_turb*np.log(yPlus[:]) + alphaHat_turb
#lamlawU = np.zeros(len(U))
#lamlawU[:] = 	yPlus			#betaHat_lam*yPlus[:] + alphaHat_lam
##
#
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.xticks(fontsize=25)
#plt.yticks(fontsize=25)
#plt.semilogx(yPlus,UPlus,marker = 'o',linestyle = ' ',color='k')
#plt.semilogx(yPlus,lamlawU,marker = ' ',linestyle = '-',color='k',linewidth='1.5')
#plt.semilogx(yPlus,loglawU,marker = ' ',linestyle = '-',color='k',linewidth='1.5')
#plt.axis([0,np.max(yPlus),0,np.max(UPlus)+0.05*np.max(UPlus)])
#plt.xlabel(r'$\mathbf{z^+}$',fontsize=30)
#plt.ylabel(r'$\mathbf{U^+}$',fontsize=30)
#plt.show()
#plt.savefig('../Data/processedData/figures/tempProfile.png')
#plt.close()










##	Plot data
#plt.plot(data.z,data.UxMean,marker='o',linestyle=' ')
#plt.show()
#plt.plot(data.z,data.UyMean,marker='o',linestyle=' ')
#plt.show()
#plt.plot(data.z,data.uxRMS,marker='o',linestyle=' ')
#plt.show()
#plt.plot(data.z,data.uyRMS,marker='o',linestyle=' ')
#plt.show()
#plt.plot(data.z,data.uv,marker='o',linestyle=' ')
#plt.show()
########################################################################################



