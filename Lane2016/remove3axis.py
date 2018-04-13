import pyfits

templatehdulist = pyfits.open('OrionA_850_auto_mos_clip.fits')
templatedata = templatehdulist[0].data[0,:,:]
templatehdulist[0].header['NAXIS'] = 2
del templatehdulist[0].header['*3']
pyfits.writeto('nofreq_OrionA_850_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

