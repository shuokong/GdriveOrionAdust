import pyfits
import numpy as np

templatehdulist = pyfits.open('OrionA_850_auto_mos_clip.fits')
#templatehdulist = pyfits.open('OrionA_850_auto_mos_clip_smooth60.fits')
#templatedata = templatehdulist[0].data[0,:,:]*26.8 # pixels per 14.6" beam
templatedata = templatehdulist[0].data[0,:,:]
#templatedata[np.isnan(templatedata)] = 0 # this is necessary for miriad convol, since convol does not handle well nan
templatehdulist[0].header['NAXIS'] = 2
#templatehdulist[0].header['BMAJ'] = 60./3600.
#templatehdulist[0].header['BMIN'] = 60./3600.
del templatehdulist[0].header['COMMENT']
del templatehdulist[0].header['']
del templatehdulist[0].header['HISTORY']
templatehdulist[0].header['OBJECT'] = 'OrionA'
for kk in templatehdulist[0].header.keys():
    if kk not in ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','OBJECT','BUNIT','DATE','TELESCOP','DATE-OBS','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2','CDELT1','CDELT2','MJD-OBS','RADESYS','EQUINOX','BACKEND']:
        del templatehdulist[0].header[kk]
del templatehdulist[0].header['LBOUND*']
#templatehdulist[0].header['BMAJ'] = 14.6/3600.
#templatehdulist[0].header['BMIN'] = 14.6/3600.
#templatehdulist[0].header['BPA'] = 0
#templatehdulist[0].header['BUNIT'] = 'JY/BEAM'
pyfits.writeto('nofreq_OrionA_850_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)
#pyfits.writeto('nofreq_OrionA_850_auto_mos_clip_smooth60.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

#templatehdulist = pyfits.open('OrionA_450_auto_mos_clip.fits')
#templatedata = templatehdulist[0].data[0,:,:]
#templatehdulist[0].header['NAXIS'] = 2
#del templatehdulist[0].header['*3']
#pyfits.writeto('nofreq_OrionA_450_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)
