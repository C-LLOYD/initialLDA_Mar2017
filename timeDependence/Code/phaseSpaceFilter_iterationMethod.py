##		THIS FUNCTION DOES NOT YET WORK - THE METHOD SEEMS
##		TO WORK BEST AFTER ONLY A SINGLE ITERATION SO UPDATE
##		THIS LATER IF NECESSARY
##
##Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv
##
#
##	Define the filtering function:
##	Input: velocity and method
##	Output: index of spikes after several loops.
def spikeLocator(U,method,writeString):
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
	writePath = writeString+'x_'+data.NXYZ[1]+'_z_'+data.NXYZ[3]+'_'+VariableName+'_phaseSpacePlot.png'
#
##
	fig= mpl.figure()
	ax = fig.gca(projection='3d')
	ax.plot_wireframe(X,Y,Z)
	ax.scatter(U[(rho<rho_e)], dU[(rho<rho_e)],d2U[(rho<rho_e)], zdir='z',color='k')
	ax.scatter(U[(rho>rho_e)], dU[(rho>rho_e)],d2U[(rho>rho_e)], zdir='z',color='r')
	ax.set_xlabel('U')
	ax.set_ylabel('dU')
	ax.set_zlabel('d2U')
	mpl.show()	
	return newSpikes


#
##	Define function
def PSF(data,method,writeString):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix();
	Uy = data.Uy.as_matrix();
	t = data.timeStamp.as_matrix()
#
##	Initialise filtered velocity fields
	UxFil = Ux
	UyFil = Uy
	tFil = t
	Nloops = 1
	notSpikes = UxFil == UxFil
	N_newSpikes = 1
	N_Spikes = 0
	while N_newSpikes > 0 and Nloops < 2:
		print("Loop number: " + str(Nloops))
		newSpikes = spikeLocator(UxFil,method,writeString)
		N_newSpikes = sum(newSpikes)
		N_Spikes = N_Spikes + N_newSpikes
		print("Number of new spikes: "+str(N_newSpikes))
		print("Total number of spikes: "+str(N_Spikes)+"= "+
			str(N_Spikes*float(100)/len(Ux))+"%")
		notSpikes = ~newSpikes
		UxFil = UxFil[notSpikes]
		tFil  = tFil[notSpikes]
		Nloops = Nloops + 1
#
##

	writePath = writeString+'x_'+data.NXYZ[1]+'_z_'+data.NXYZ[3]+'_'+variableName+'_filteredTimeSeries.png'
#
##
	mpl.subplot(2, 1, 1)
	mpl.plot(t,Ux)
	mpl.subplot(2, 1, 2)
	mpl.plot(tFil,UxFil)
	mpl.show()
#
##	SO:	Function currently works BUT it is harsh when removing 
##		spikes after the first iteration ...
##		For now we loop once, just to remove the largest spikes.
#
##	Now plot the result ..
		
#
##	Loop through both velocity components
#
#
##	Compute mean/median and offset U
##	According to Wahl (reference) taking median and MADS results
##	in more stability but current tests suggest mean is more
##	appropriate for us.


#
##	Plot the two data sets as time series:
###	Append the dataframe with the filtered data
###	Replace the spikes with Nan - this will allow the filtered data to be plotted against the
##	same time vector.
#		U_fil = U
#		index = rhoe>rho_e
#		U_fil[index] = np.nan
#		This needs fixing,
#		then we need to start looping,
#		Then we need to save the final files ...
#		!! we will need to update the spikes each loop too, otherwise the original spikes
#		will not be included!!
#

data = pd.read_pickle('../Data/processedData/dataFrames/8hz_x_400.0000_z_40.0000_data_weighted.pkl')
U_raw = data.Ux.as_matrix()
#spikeLocator(U_raw,'mean')

PSF(data,'median')
#
#
#
#
