import sys
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('planck_herschel_plane1.fits'))[0]
#from reproject import reproject_exact
#array, footprint = reproject_exact(hdu2, hdu1.header)
#fits.writeto('lombardi_on_stutz_header.fits', array, hdu1.header, clobber=True)

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('emap_Orion_A_bw2.0.fits'))[1]
#from reproject import reproject_exact
#array, footprint = reproject_exact(hdu2, hdu1.header)
#fits.writeto('stefan_on_stutz_header.fits', array, hdu1.header, clobber=True)

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#from reproject import reproject_exact
#array, footprint = reproject_exact(hdu2, hdu1.header)
#fits.writeto('Lane_on_Stutz_header.fits', array, hdu1.header, clobber=True)

hdu1 = fits.open(get_pkg_data_filename('emap_Orion_A_bw1.0.fits'))[1]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
hdu2 = fits.open(get_pkg_data_filename('Lane2016/convol60_nofreq_OrionA_850_auto_mos_clip.fits'))[0]
from reproject import reproject_interp
array, footprint = reproject_interp(hdu2, hdu1.header)
fits.writeto('Lane_on_Stefan_header.fits', array, hdu1.header, clobber=True)
