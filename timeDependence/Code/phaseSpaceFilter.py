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
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv
##
#
data = pd.read_pickle('../Data/processedData/dataFrames/8hz_x_400.0000_z_40.0000_data_weighted.pkl')
t = U.timeStamp.as_matrix()
##	Define function
def PSF(data):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp.as_matrix()
	U = [Ux,Uy]
#
##	Compute median and offset U
U = U - np.mean(U)
##	Compute derivatives:
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
##	Compute MADS and estimate maximum expected values#
temp = 1.483*(np.median(np.absolute(U-np.median(U))))
lambdaU = -np.sqrt(2)*erfinv(1-1/(2*len(U)))*temp
lambdaU = np.sqrt(2*np.log(len(U)))*np.std(U)
#
temp = 1.483*(np.median(np.absolute(dU-np.median(dU))))
lambdadU = -np.sqrt(2)*erfinv(1-1/(2*len(dU)))*temp
lambdadU = np.sqrt(2*np.log(len(dU)))*np.std(dU)
#
temp = 1.483*(np.median(np.absolute(d2U-np.median(d2U))))
lambdad2U = -np.sqrt(2)*erfinv(1-1/(2*len(d2U)))*temp
lambdad2U = np.sqrt(2*np.log(len(d2U)))*np.std(d2U)
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
rho_e1 = (
	np.divide((np.sin(phi)*np.cos(theta)*np.cos(alpha) + np.cos(phi)*np.sin(alpha))**2,a**2)
	+
	np.divide((np.sin(phi)*np.cos(theta)*np.sin(alpha) + np.cos(phi)*np.cos(alpha))**2,b**2)
	+
	np.divide((np.sin(phi)*np.sin(theta))**2,c**2)
		)**(-0.5)


rho_e = (
	np.divide((np.sin(phi)*np.cos(theta))**2,a**2)
	+
	np.divide((np.sin(phi)*np.sin(theta))**2,b**2)
	+
	np.divide((np.cos(phi))**2,c**2)
		)**(-0.5)
##	Compare lengths rho and rho_e and identify spikes
print(len(U[(rho>rho_e1)]))
print(len(U[(rho>rho_e)]))

spikesU = U[(rho>rho_e)]
spikesdU = dU[(rho>rho_e)]
spikesd2U = d2U[(rho>rho_e)]
#
##	Convert back to cart. to plot ellipse
N=len(theta)
#N=27

#theta = np.linspace(0,2*np.pi,N)
#phi = np.linspace(0,np.pi,N)

X = a*(np.outer(np.sin(phi),np.cos(theta))*np.cos(alpha)+np.outer(np.cos(phi),np.ones(N))*np.sin(alpha))
Y = b*(np.outer(np.sin(phi),np.cos(theta))*np.sin(alpha)-np.outer(np.cos(phi),np.ones(N))*np.cos(alpha))
Z = c*(np.outer(np.sin(phi),np.sin(theta)))


fig= mpl.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(X,Y,Z)
#ax.scatter(U[(rho<rho_e)], dU[(rho<rho_e)],d2U[(rho<rho_e)], zdir='z',color='k')
ax.scatter(U[(rho>rho_e)], dU[(rho>rho_e)],d2U[(rho>rho_e)], zdir='z',color='r')
ax.set_xlabel('U')
ax.set_ylabel('dU')
ax.set_zlabel('d2U')
mpl.show()

print(len(spikesU))
print(np.std(U))
print(1.483*(np.median(np.absolute(U-np.median(U)))))
print(np.std(U[(rho<rho_e)]))
print(1.483*(np.median(np.absolute(U[(rho<rho_e)]-np.median(U[(rho<rho_e)])))))

#fig= mpl.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(U, dU,d2U, zdir='z')
#ax.scatter(spikesU, spikesdU,spikesd2U, zdir='z',color='r')
#mpl.show()
#
#fig= mpl.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(U[(rho<rho_e)], dU[(rho<rho_e)],d2U[(rho<rho_e)], zdir='z')
#ax.plot_wireframe(X,Y,Z)
#mpl.show()
#
#fig= mpl.figure()
#ax = fig.add_subplot(111, projection='3d')
##ax = fig.gca(projection='3d')
#ax.scatter(X,Y,Z, zdir='z')
#mpl.show()
#
#fig= mpl.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(U[(rho<rho_e)], dU[(rho<rho_e)],d2U[(rho<rho_e)], zdir='z')
#mpl.show()
#
#
#
mpl.subplot(2, 1, 1)
mpl.plot(t,U)
mpl.scatter(t[(rho>rho_e)],U[(rho>rho_e)],color='r',s=20)
mpl.axis([0,300,-0.15,0.15])
mpl.subplot(2, 1, 2)
mpl.plot(t[(rho<rho_e)],U[(rho<rho_e)])
mpl.axis([0,300,-0.15,0.15])
mpl.show()

