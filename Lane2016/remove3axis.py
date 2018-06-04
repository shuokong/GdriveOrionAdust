import pyfits

templatehdulist = pyfits.open('OrionA_850_auto_mos_clip.fits')
templatedata = templatehdulist[0].data[0,:,:]*3.1415926*(14.6/2.)**2/3.**2
templatehdulist[0].header['NAXIS'] = 2
templatehdulist[0].header['BMAJ'] = 14.6/3600.
templatehdulist[0].header['BMIN'] = 14.6/3600.
templatehdulist[0].header['BPA'] = 0
templatehdulist[0].header['BUNIT'] = 'JY/BEAM'
del templatehdulist[0].header['COMMENT']
del templatehdulist[0].header['']
del templatehdulist[0].header['HISTORY']
templatehdulist[0].header['OBJECT'] = 'OrionA'
for kk in templatehdulist[0].header.keys():
    if kk not in ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','OBJECT','BUNIT','DATE','TELESCOP','DATE-OBS','CRPIX1','CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2','CDELT1','CDELT2','MJD-OBS','RADESYS','EQUINOX','BACKEND','BMAJ','BMIN','BPA']:
        del templatehdulist[0].header[kk]
#del templatehdulist[0].header['LBOUND*']
pyfits.writeto('nofreq_OrionA_850_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

#templatehdulist = pyfits.open('OrionA_450_auto_mos_clip.fits')
#templatedata = templatehdulist[0].data[0,:,:]
#templatehdulist[0].header['NAXIS'] = 2
#del templatehdulist[0].header['*3']
#pyfits.writeto('nofreq_OrionA_450_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)
