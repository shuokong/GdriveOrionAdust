import sys
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('planck_herschel_plane1.fits'))[0]
#from reproject import reproject_exact
#array, footprint = reproject_exact(hdu2, hdu1.header)
#fits.writeto('lombardi_on_stutz_header.fits', array, hdu1.header, clobber=True)

hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('emap_Orion_A_bw2.0.fits'))[1]
from reproject import reproject_exact
array, footprint = reproject_exact(hdu2, hdu1.header)
fits.writeto('stefan_on_stutz_header.fits', array, hdu1.header, clobber=True)


