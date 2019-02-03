# from https://pyregion.readthedocs.io/en/latest/getting_started.html

import pyregion
from astropy.io import fits
import numpy as np

region_name = "OrionKLellipse.reg"
r = pyregion.open(region_name)

hdu1 = fits.open('Lane_on_Stefan_header_CASA.fits')[0]
mymask = r.get_mask(hdu=hdu1)
hdu1.data[mymask] = np.nan
fits.writeto('OrionKLellipse_Lane_on_Stefan_header_CASA.fits', hdu1.data, hdu1.header, clobber=True)

hdu1 = fits.open('mask_emap_Orion_A_bw1.0.fits')[0]
mymask = r.get_mask(hdu=hdu1)
hdu1.data[mymask] = np.nan
fits.writeto('OrionKLellipse_mask_emap_Orion_A_bw1.0.fits', hdu1.data, hdu1.header, clobber=True)

hdu1 = fits.open('mask_emap_Orion_A_bw1.0.fits')[0]
hdu1.data[~np.isnan(hdu1.data)] = 1
hdu1.data[np.isnan(hdu1.data)] = 0
fits.writeto('bool_mask_emap_Orion_A_bw1.0.fits', hdu1.data, hdu1.header, clobber=True)
