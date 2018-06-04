import pyfits
import sys

templatehdulist = pyfits.open(sys.argv[1])
print templatehdulist.info()
neg = templatehdulist[0].data < 0
templatehdulist[0].data[neg] = 0
templatedata = templatehdulist[0].data[int(sys.argv[2]),:,:]
templatehdulist[0].header['NAXIS'] = 2
#del templatehdulist[0].header['CRPIX3']
#del templatehdulist[0].header['CDELT3']
#del templatehdulist[0].header['CRVAL3']
#del templatehdulist[0].header['CTYPE3']
del templatehdulist[0].header['NAXIS3']
del templatehdulist[0].header['PLANE2']
del templatehdulist[0].header['PLANE3']
del templatehdulist[0].header['PLANE4']
del templatehdulist[0].header['PLANE5']
del templatehdulist[0].header['PLANE6']
del templatehdulist[0].header['PLANE7']
pyfits.writeto(sys.argv[3],templatedata,templatehdulist[0].header,output_verify='exception',clobber=True,checksum=False)

