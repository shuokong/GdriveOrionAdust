import sys
import numpy as np
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
#hdu2 = fits.open(get_pkg_data_filename('emap_Orion_A_bw1.0.fits'))[1]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu1, hdu2.header)
#fits.writeto('stutz_on_stefan_header.fits', array, hdu2.header, clobber=True)

hdu2 = fits.PrimaryHDU()
hdu2.header['BITPIX']  = -32 
hdu2.header['NAXIS']   = 2
hdu2.header['NAXIS1']  = 560
hdu2.header['NAXIS2']  = 292
hdu2.header['CTYPE1']  = 'GLON-TAN'                                                            
hdu2.header['CTYPE2']  = 'GLAT-TAN'                                                            
hdu2.header['CRVAL1']  = 210.09792
hdu2.header['CRVAL2']  = -19.912362
hdu2.header['CUNIT1']  = 'deg'                                                                 
hdu2.header['CUNIT2']  = 'deg'                                                                 
hdu2.header['CDELT1']  = -0.002948047285127
hdu2.header['CDELT2']  = 0.002948047285127
hdu2.header['COORDSYS'] = 'gal'                                                                 
hdu2.header['CRPIX1']  = 1
hdu2.header['CRPIX2']  = 1
hdu1 = fits.open('../AncillaryData/GBT/OrionA_Tkin_DR1_rebase3_flag.fits')[0]
from reproject import reproject_interp
array, footprint = reproject_interp(hdu1, hdu2.header)
fits.writeto('../AncillaryData/GBT/ga_OrionA_Tkin_DR1_rebase3_flag.fits', array, hdu2.header, clobber=True)

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/smooth18_nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu2, hdu1.header)
#del hdu1.header['COMMENT']
#del hdu1.header['HISTORY']
#fits.writeto('Lane_on_Stutz_header.fits', array, hdu1.header, clobber=True)

#hdu1 = fits.open(get_pkg_data_filename('emap_Orion_A_bw1.0.fits'))[1]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip_smooth60.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/convol60_OrionA_850_auto_mos_clip.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lanecores_nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu2, hdu1.header)
#fits.writeto('Lane_on_Stefan_header.fits', array, hdu1.header, clobber=True)
#fits.writeto('Lane_on_Stefan_header_Starlink.fits', array, hdu1.header, clobber=True)
#fits.writeto('Lane_on_Stefan_header_CASA.fits', array, hdu1.header, clobber=True)
#hdu1.data[np.isnan(array)] = np.nan
#fits.writeto('mask_emap_Orion_A_bw1.0.fits', hdu1.data, hdu1.header, clobber=True)
#fits.writeto('Lanecores_on_Stefan_header.fits', array, hdu1.header, clobber=True) # not good

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA-all_conv500_temp.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('OrionKLellipse_Lane_on_Stefan_header_CASA.fits'))[0]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu1, hdu2.header)
#array[np.isnan(hdu2.data)] = np.nan
#fits.writeto('dusttemp_Stutz_on_OrionKLellipse_header.fits', array, hdu2.header, clobber=True)

#hdu1 = fits.open(get_pkg_data_filename('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu1, hdu2.header)
#array[np.isnan(hdu2.data)] = np.nan
#fits.writeto('herschelAmelia/Stutz_on_Lane_header.fits', array, hdu2.header, clobber=True)
#dist = 400. # pc
#pixelscale = 3.*dist*1.5e13 # cm
#print 'np.nansum(array)',np.nansum(array)
#print 'total Herschel mass in Lane region',np.nansum(array)/4.27e23*pixelscale**2/2.e33,'Msun'

#hdu1 = fits.open(get_pkg_data_filename('emap_Orion_A_bw1.0.fits'))[1]
#hdu2 = fits.open(get_pkg_data_filename('Lane2016/nofreq_OrionA_850_auto_mos_clip.fits'))[0]
#from reproject import reproject_interp
#array, footprint = reproject_interp(hdu1, hdu2.header)
#fits.writeto('Stefan_on_Lane_header.fits', array, hdu2.header, clobber=True) 

