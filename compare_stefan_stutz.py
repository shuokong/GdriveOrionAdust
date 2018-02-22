import sys
import os
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
font = {'weight' : 'normal','size':14,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('stefan_on_stutz_header.fits'))[0]
stefan=hdu2.data
stutz=hdu1.data

p=plt.figure(figsize=(7,6))
plt.subplots_adjust(top=0.97,bottom=0.12,left=0.12,right=0.96)
ax=p.add_subplot(111)
plt.plot(stutz,stefan,'k.')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'stutz $N_H$')
plt.ylabel(r'stefan $A_K$')
#plt.xlim(0,1e24)

os.system('rm stefan_stutz.png')
plt.savefig('stefan_stutz.png')
os.system('open stefan_stutz.png')

selection = (stutz>1e23)
stefan[selection] = 1
stefan[np.invert(selection)] = 0
hdu2.writeto('selection.fits',overwrite=True)


