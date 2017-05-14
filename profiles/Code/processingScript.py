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
import matplotlib.pyplot as mpl
import pandas as pd
import glob
from processingFunctions import txtToDataFrame
from processingFunctions import timeAverage
from processingFunctions import plotter
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
rawPath = 	["../Data/rawData/4hz/300mm/*/*.txt",
			"../Data/rawData/4hz/400mm/*/*.txt",
			"../Data/rawData/8hz/300mm/*/*.txt",
			"../Data/rawData/8hz/400mm/*/*.txt"]
#
errorPath = 	["../Data/processedData/dataFrames/4hz_errors.pkl",
			"../Data/processedData/dataFrames/4hz_errors.pkl",
			"../Data/processedData/dataFrames/8hz_errors.pkl",
			"../Data/processedData/dataFrames/8hz_errors.pkl"]
#
writePath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl"]
#
##
dataPath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl"]
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
#
##	Rearrange and save dataFrame - if wanted
##	Rearrange does not work as intended : 1,10,2, ...
		data = data.sort_values(by='z',ascending=1)
#
##	Now add columns for the errors, reading from the errorPath
#
#
##
		data = data.set_index("z")	
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
		data = data.reset_index(level=0)
#
#	print(data)
		data.to_pickle(writePath[j])


#mpl.semilogx(data.z,data.UxMean,linestyle = 'None', marker = '.')
#mpl.errorbar(data.z,data.UxMean,yerr=data.error_UxMean, linestyle = "None")
#mpl.show()
#mpl.semilogx(data.z,data.UyMean,linestyle = 'None', marker = '.')
#mpl.errorbar(data.z,data.UyMean,yerr=data.error_UyMean, linestyle = "None")
#mpl.show()
#mpl.semilogx(data.z,data.uxRMS,linestyle = 'None', marker = '.')
#mpl.errorbar(data.z,data.uxRMS,yerr=data.error_uxRMS, linestyle = "None")
#mpl.show()
#mpl.semilogx(data.z,data.uyRMS,linestyle = 'None', marker = '.')
#mpl.errorbar(data.z,data.uyRMS,yerr=data.error_uyRMS, linestyle = "None")
#mpl.show()
#mpl.semilogx(data.z,data.uv,linestyle = 'None', marker = '.')
#mpl.errorbar(data.z,data.uv,yerr=data.error_uv, linestyle = "None")
#mpl.show()
#
#
###############	TEMP LOCATION OF LAW OF THE WALL FUNCTION - MOVE THIS LATER ################
####
####	1. 	define new variables
##	Initialise variables using the full data set and iterate to find log-law region by removing extremes of z
def linearReg(X,Y):
	beta = np.sum( (X-np.mean(X)) * (Y-np.mean(Y)) ) / np.sum( (Y - np.mean(Y))**2 ) 
	alpha = np.mean(X) - beta*np.mean(Y)
	return [alpha,beta]

data = pd.read_pickle(dataPath[0])
Ulab = data.UxMean.as_matrix()
Ylab = data.z.as_matrix()/1000
nu = 1.1*10**(-6)
kappa = 0.41
Cplus = 5.2
#
##
Ylam = Ylab[0:20]
Ulam = Ylab[0:20]
counter = 0
while counter < 5:
	counter = counter + 1
#	print(counter)
#	print(Ylam,Ulam)
	[aLam,bLam] = linearReg(Ulam,Ylam)
#	Delta = aLam/bLam
	Utau1 = np.sqrt(bLam*nu)
	Delta = aLam*nu/(Utau1**2)
	Yplus = (Ylab+Delta)*Utau1/nu
#
##	Use these to estimate range for log region
	Yturb = Ylab[(30<=Yplus) & (Yplus<=200)]
	Uturb = Ulab[(30<=Yplus) & (Yplus<=200)]
#
##	Now estimate Utau based on Delta and linear regression of log law component
	[aturb,bTurb] = linearReg(Uturb,np.log(Yturb+Delta))
	Utau2 = bTurb*kappa
	Yplus = Ylab*Utau2/nu
	Ylam = Ylab[(Yplus<=5)]
#	print(len(Ylam))
	Ulam = Ulab[(Yplus<=5)]
	if len(Ylam) < 2:
		Ylam = Ylam[0:3]
		Ulam = Ulam[0:3]
#	print(Utau1,Utau2)

Yplus = (Ylab + Delta)*Utau1/nu
Uplus = Ulab/Utau1

lamlawUplus = Yplus
loglawUplus = np.log(Yplus)/kappa + Cplus + (Utau1/kappa)*np.log(Utau1/nu)

mpl.semilogx(Yplus,Uplus,linestyle=' ',marker = '.')
mpl.semilogx(Yplus,lamlawUplus)
mpl.semilogx(Yplus,loglawUplus)
mpl.axis([0,np.max(Yplus),0,np.max(Uplus)+0.05*np.max(Uplus)])
mpl.show()
mpl.close()

mpl.plot(Uplus,Yplus,linestyle=' ',marker= '.')
mpl.plot(lamlawUplus,Yplus)
mpl.plot(loglawUplus,Yplus)	
mpl.axis([0,30,0,100])
mpl.show()
mpl.close()	

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
#mpl.rc('text', usetex=True)
#mpl.rc('font', family='serif')
#mpl.xticks(fontsize=25)
#mpl.yticks(fontsize=25)
#mpl.semilogx(yPlus,UPlus,marker = 'o',linestyle = ' ',color='k')
#mpl.semilogx(yPlus,lamlawU,marker = ' ',linestyle = '-',color='k',linewidth='1.5')
#mpl.semilogx(yPlus,loglawU,marker = ' ',linestyle = '-',color='k',linewidth='1.5')
#mpl.axis([0,np.max(yPlus),0,np.max(UPlus)+0.05*np.max(UPlus)])
#mpl.xlabel(r'$\mathbf{z^+}$',fontsize=30)
#mpl.ylabel(r'$\mathbf{U^+}$',fontsize=30)
#mpl.show()
#mpl.savefig('../Data/processedData/figures/tempProfile.png')
#mpl.close()










##	Plot data
#mpl.plot(data.z,data.UxMean,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.UyMean,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uxRMS,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uyRMS,marker='o',linestyle=' ')
#mpl.show()
#mpl.plot(data.z,data.uv,marker='o',linestyle=' ')
#mpl.show()
########################################################################################



