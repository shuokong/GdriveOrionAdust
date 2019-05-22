import sys
import os
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'normal','size':30,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

#hdu1 = fits.open(get_pkg_data_filename('stutz_on_stefan_header.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('emap_Orion_A_bw1.0.fits'))[1]
#pdfname = 'stefan_stutz_1arcmin.pdf'
hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('stefan_on_stutz_header.fits'))[0]
pdfname = 'stefan_stutz_18arcsec.pdf'
stefan=hdu2.data
stutz=hdu1.data
selection = (stutz>2e19)&(stutz<2e24)&(stefan>1e-4)&(stefan<1e1)

p=plt.figure(figsize=(16,12))
plt.subplots_adjust(top=0.97,bottom=0.08,left=0.09,right=0.98)
ax=p.add_subplot(111)
plt.plot(stutz[selection],stefan[selection],'k.', rasterized=True)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$N_H~\rm (cm^{-2})$')
plt.ylabel(r'$A_K~\rm (mag)$')
#plt.xlim(0,1e24)
plt.ylim(1e-4,10)
plt.vlines(1e23,1e-4,10,'k',linestyles='dashed')
plt.plot([np.nanmin(stutz[selection]),np.nanmax(stutz[selection])],[np.nanmin(stutz[selection])/2.e21/8.93,np.nanmax(stutz[selection])/2.e21/8.93],'g--')

os.system('rm '+pdfname)
plt.savefig(pdfname)
os.system('open '+pdfname)
os.system('cp '+pdfname+' ~/GoogleDrive/imagesSFE/')

#selection = (stutz>1e23)
#stefan[selection] = 1
#stefan[np.invert(selection)] = 0
#hdu2.writeto('selection.fits',overwrite=True)


