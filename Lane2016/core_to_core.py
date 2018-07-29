import sys
import os
from scipy.optimize import curve_fit
import numpy as np
import telescope_sensitivity as ts
import pyfits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
from astropy.io import fits
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':16,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

def gaus(x,a,x0,sigma):
    return a/sigma*np.exp(-(x-x0)**2/(2*sigma**2))

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

yso = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/GAScores.txt',usecols=(6),unpack=True,dtype='string')
ysoyes = (yso == 'Y')
ysono = (yso == 'N')
#savetxtarr = np.stack((corenames,xw,yw,coremasses,corevelocities[0],coresnr[0],corevelocities[1],coresnr[1],corevelocities[2],coresnr[2],nh3velocities[0]),axis=1)
corenames,xw,yw,coremasses,corevelocities12CO,coresnr12CO,coremom0_12CO,coremom1_12CO,coremom2_12CO,corevelocities13CO,coresnr13CO,coremom0_13CO,coremom1_13CO,coremom2_13CO,corevelocitiesC18O,coresnrC18O,coremom0_C18O,coremom1_C18O,coremom2_C18O,nh3velocities,nh3evelocities,nh3tkin,coregaus_x0_13CO,coregaus_ex0_13CO,coregaus_sigma_13CO,coregaus_esigma_13CO,coregaus_x0_C18O,coregaus_ex0_C18O,coregaus_sigma_C18O,coregaus_esigma_C18O = np.loadtxt('convol32_Kirkcores_velocities.txt',unpack=True)
#print np.stack((corevelocities13CO.T,corevelocitiesC18O.T),axis=1)
coreysoyes = corenames[ysoyes]
coreysono = corenames[ysono]
print len(coreysoyes),len(coreysono)
#sys.exit()
worldcoord2 = np.stack((xw,yw),axis=1) # 

nhhdu = fits.open('../herschelAmelia/carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
maskhdu = fits.open('../herschelAmelia/mask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
north_maskhdu = fits.open('../herschelAmelia/northmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
south_maskhdu = fits.open('../herschelAmelia/southmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]

w = wcs.WCS(nhhdu.header)
pixcoord = w.all_world2pix(worldcoord2,1) # FITS standard uses 1
xx = pixcoord[:,0]
yy = pixcoord[:,1]

topcluster_velocity = [nh3velocities[nn] for nn,ii in enumerate(corenames) if north_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1]
print 'topcluster_velocity[0,1,2]',topcluster_velocity[:3]
botcluster_velocity = [nh3velocities[nn] for nn,ii in enumerate(corenames) if south_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1]
print 'botcluster_velocity[0,1,2]',botcluster_velocity[:3]
totcluster_velocity = [nh3velocities[nn] for nn,ii in enumerate(corenames) if maskhdu.data[int(yy[nn]),int(xx[nn])] == 1]

print 'np.nanstd(topcluster_velocity)',np.nanstd(topcluster_velocity) 
print 'np.nanstd(botcluster_velocity)',np.nanstd(botcluster_velocity) 
print 'np.nanstd(totcluster_velocity)',np.nanstd(totcluster_velocity) 

dist = 400. # pc
pixelscale = 6.*dist*1.5e13 # cm

totmass = np.nansum(nhhdu.data*maskhdu.data)/4.27e23*pixelscale**2
north_mass = np.nansum(nhhdu.data*north_maskhdu.data)/4.27e23*pixelscale**2
south_mass = np.nansum(nhhdu.data*south_maskhdu.data)/4.27e23*pixelscale**2
print 'totmass in Msun',totmass/2.e33

filength = 34.*3.*60.*dist*1.5e13 # 34 segments in ISF each 3' like data paper PV cut
print 'tot mass per unit length (Msun/pc)',(totmass/2.e33)/(filength/3.08e18)
north_filength = 9.*3.*60.*dist*1.5e13 # 
print 'north mass per unit length Msun/pc',(north_mass/2.e33)/(north_filength/3.08e18)
south_filength = 20.*3.*60.*dist*1.5e13 # 
print 'south mass per unit length Msun/pc',(south_mass/2.e33)/(south_filength/3.08e18)

sigma_virial = (6.67e-8*totmass/filength/2.)**0.5/1.e5
print 'sigma_virial',sigma_virial
north_sigma_virial = (6.67e-8*north_mass/north_filength/2.)**0.5/1.e5
print 'north_sigma_virial',north_sigma_virial
south_sigma_virial = (6.67e-8*south_mass/south_filength/2.)**0.5/1.e5
print 'south_sigma_virial',south_sigma_virial

