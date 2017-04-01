#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
import numpy as np
import re
import matplotlib.pyplot as mpl
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv

#	Import data:
U = pd.read_pickle('../Data/processedData/dataFrames/4hz_x_400.0000_z_1.0000_data_weighted.pkl')
U = U.Uy.as_matrix()
#
##	Compute median and offset U
U = U - np.median(U)
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
#
temp = 1.483*(np.median(np.absolute(dU-np.median(dU))))
lambdadU = -np.sqrt(2)*erfinv(1-1/(2*len(dU)))*temp
#
temp = 1.483*(np.median(np.absolute(d2U-np.median(d2U))))
lambdad2U = -np.sqrt(2)*erfinv(1-1/(2*len(d2U)))*temp
#
##	Compute rotation angle of the principle axis of d2U vs U using cross correlation
alpha = np.arctan(sum(U*d2U)/sum(U**2))
##	Now compute constants a, b and c for the ellipsoid definition
b = np.sqrt(np.divide(
	lambdad2U**2 - (lambdaU**2)*(np.sin(alpha))**2,
	(np.cos(alpha))**2 - ((np.tan(alpha))**2)*((np.sin(alpha))**2)
		))
#
a = np.sqrt(lambdaU**2 - (b**2)*(np.cos(theta))**2)
#
c = lampdadU


#fig= mpl.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(U, dU,d2U, zdir='z')
#mpl.show()

