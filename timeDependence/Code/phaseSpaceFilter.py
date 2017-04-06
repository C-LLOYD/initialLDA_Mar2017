#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv
##
#
##	Define the filtering function:
##	Input: velocity and method
##	Output: index of spikes after several loops.
def spikeLocator(U,data,method,writePaths_figures,VariableName):
	if 	method == 'mean':
		U = U - np.mean(U)
	elif 	method == 'median':
		U = U - np.median(U)
	else:
		print('No averaging technique given')
#
##	Compute derivatives:
##	Central differencing used for central elements and gradient 
##	assumed constant over first and last elements.
	dU = np.zeros(len(U))
	dU[:] = np.NAN
	dU[1:-1] = (U[2:]-U[0:-2])/2
	dU[0] = dU[1]; dU[-1] = dU[-2];
#
	d2U = np.zeros(len(U))
	d2U[:] = np.NAN
	d2U[1:-1] = (dU[2:]-dU[0:-2])/2
	d2U[0] = d2U[1]; d2U[-1] = d2U[-2];
#
##	Use these to calculate maximum expected values:
##	TWO METHODS:
##		STD + MEAN or MAD + MEDIAN
	if 	method == 'mean':
		lambdaU = np.sqrt(2*np.log(len(U)))*np.std(U)
		lambdadU = np.sqrt(2*np.log(len(dU)))*np.std(dU)
		lambdad2U = np.sqrt(2*np.log(len(d2U)))*np.std(d2U)
#
	elif 	method == 'median':
		temp = 1.483*(np.median(np.absolute(U-np.median(U))))
		lambdaU = -np.sqrt(2)*erfinv(1-(2*len(U))**(-1))*temp
		temp = 1.483*(np.median(np.absolute(dU-np.median(dU))))
		lambdadU = -np.sqrt(2)*erfinv(1-(2*len(dU))**(-1))*temp
		temp = 1.483*(np.median(np.absolute(d2U-np.median(d2U))))
		lambdad2U = -np.sqrt(2)*erfinv(1-(2*len(d2U))**(-1))*temp
#
	else:
		print('No averaging technique given')
#
##	Compute rotation angle of the principle axis of d2U vs U using cross correlation
	alpha = np.arctan2(sum(U*d2U),sum(U**2))
##	Now compute constants a, b and c for the ellipsoid definition
	b = np.sqrt(np.divide(
		lambdad2U**2 - (lambdaU**2)*(np.sin(alpha))**2,
		(np.cos(alpha))**2 - ((np.tan(alpha))**2)*((np.sin(alpha))**2)
		))
#
	a = np.sqrt(lambdaU**2 - (b**2)*(np.cos(alpha))**2)
#
	c = lambdadU
#
##	Convert to sph coord
	theta = np.arctan2(dU,U)
	phi = np.arctan2(np.sqrt(U**2 + dU**2),d2U)
	rho = np.sqrt(U**2 + dU**2 + d2U**2)
#
##	Compute rho_e based on theta and phi
##	This gives the location of the ellipsoid surface for each angle 
##	which can be compared against the actual distance rho.
##	This definition of rho_e is based on the definition of an ellipsoid (wolfram alpha / textbooks)
	rho_e = (
		np.divide((np.sin(phi)*np.cos(theta))**2,a**2)
		+
		np.divide((np.sin(phi)*np.sin(theta))**2,b**2)
		+
		np.divide((np.cos(phi))**2,c**2)
		)**(-0.5)
#
##	Compare lengths rho and rho_e and identify spikes
##	All occurances of rho>rho_e give us the poor data, thus the function could be exited here.
##	Else, we can plot the data in phase space to provide examples of the filtering function
##	working.
##
##	Currently we have no method to suggest the filter isn't working properly (i.e if the 
##	filter suggests 90% of data is poor)
##	Add this functionality in at a later date.
	newSpikes = rho>rho_e	
##	Now plot the data: Two/ three plots
##	1.	phase space with all data and ellipsoid
##	2.	two sub plots of time series data, with and without spikes removed
#
##	Plot phase space:
##	Convert back to cart. to plot ellipse
	N=27
#
	theta = np.linspace(0,2*np.pi,N)
	phi = np.linspace(0,np.pi,N)
#
##	Use definition from Wahl (reference)
	X = a*(np.outer(np.sin(phi),np.cos(theta))*np.cos(alpha)+np.outer(np.cos(phi),np.ones(N))*np.sin(alpha))
	Y = b*(np.outer(np.sin(phi),np.cos(theta))*np.sin(alpha)-np.outer(np.cos(phi),np.ones(N))*np.cos(alpha))
	Z = c*(np.outer(np.sin(phi),np.sin(theta)))
#
##
	writePath = writePaths_figures+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_'+VariableName+'_phaseSpacePlot.png'
#
##
	fig= mpl.figure()
	mpl.rc('text', usetex=True)
	ax = fig.gca(projection='3d')
	ax.plot_wireframe(X,Y,Z)
	ax.scatter(U[(rho<rho_e)], dU[(rho<rho_e)],d2U[(rho<rho_e)], zdir='z',color='k')
	ax.scatter(U[(rho>rho_e)], dU[(rho>rho_e)],d2U[(rho>rho_e)], zdir='z',color='r')
	ax.set_xlabel(r'$\mathbf{U}$ (m/s)',fontsize=20)
	ax.xaxis.labelpad=15
	ax.set_ylabel(r'$\mathbf{d U}$ (m/s)',fontsize=20)
	ax.yaxis.labelpad=15
	ax.set_zlabel(r'$\mathbf{d^2 U}$ (m/s)',fontsize=20)
	ax.zaxis.labelpad=15
##	We could save the plot but currently this saves on every loop: Not useful!
##	Have it save on only the first loop or not at all.
	mpl.savefig(writePath)
#	mpl.show()
	mpl.close()
#
	return newSpikes

#
##	Define function
def phaseSpaceFilter(data,method,writePaths_figures,writePath_dataFrames):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp.as_matrix()
	s = data.sampleNumber.as_matrix()
	resT = data.resTime.as_matrix()
#
##	Initialise filtered velocity fields
	XSpikes = spikeLocator(Ux,data,method,writePaths_figures,'Ux')
	YSpikes = spikeLocator(Uy,data,method,writePaths_figures,'Uy')
	Spikes = XSpikes + YSpikes
	N_Spikes = sum(Spikes)
	print("Number of Ux spikes: "+str(sum(XSpikes)))
	print("Number of Uy spikes: "+str(sum(YSpikes)))
	print("Total number of spikes: "+str(N_Spikes)+" = "+
		str(N_Spikes*float(100)/len(Ux))+"%")
	Ux = Ux[~Spikes]
	Uy = Uy[~Spikes]
	t  = t[~Spikes]
	resT  = resT[~Spikes]
	s  = s[~Spikes]
#
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	Ux = pd.Series(Ux)
	Uy = pd.Series(Uy)
	t = pd.Series(t)
	s = pd.Series(s)
	resT = pd.Series(resT)
#
	data['Ux']= Ux
	data['Uy']= Uy
	data['timeStamp']= t
	data['sampleNumber']=s
	data['resTime']=resT
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
#	data.to_pickle(writePath_dataFrames+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_data_filtered.pkl')
	return data;
##
##
##
