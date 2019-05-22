import pyfits

templatehdulist = pyfits.open('OrionA_850_auto_mos_clip.fits')
templatedata = templatehdulist[0].data*1.133*14.6**2/3.**2
del templatehdulist[0].header['HISTORY']
templatehdulist[0].header['BMAJ'] = 14.6/3600.
templatehdulist[0].header['BMIN'] = 14.6/3600.
templatehdulist[0].header['BPA'] = 0
templatehdulist[0].header['BUNIT'] = 'JY/BEAM'
pyfits.writeto('beam_OrionA_850_auto_mos_clip.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

