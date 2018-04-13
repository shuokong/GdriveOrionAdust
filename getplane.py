import pyfits

templatehdulist = pyfits.open('planck_herschel.fits')
neg = templatehdulist[0].data < 0
templatehdulist[0].data[neg] = 0
#templatedata = templatehdulist[0].data[0,:,:]
#templatedata = templatehdulist[0].data[2,:,:]
templatedata = templatehdulist[0].data[3,:,:]
templatehdulist[0].header['NAXIS'] = 2
del templatehdulist[0].header['NAXIS3']
del templatehdulist[0].header['PLANE1']
del templatehdulist[0].header['PLANE2']
del templatehdulist[0].header['PLANE3']
#del templatehdulist[0].header['PLANE4']
del templatehdulist[0].header['PLANE5']
del templatehdulist[0].header['PLANE6']
del templatehdulist[0].header['PLANE7']
#pyfits.writeto('planck_herschel_plane1.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)
#pyfits.writeto('lombardi_planck_herschel_plane3_colorT.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)
pyfits.writeto('lombardi_planck_herschel_plane4_colorTerror.fits',templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

