# from https://pyregion.readthedocs.io/en/latest/getting_started.html

import pyregion
from astropy.io import fits

region_name = "totfil.reg"
#region_name = "northfil.reg"
#region_name = "southfil.reg"
r = pyregion.open(region_name)

hdu1 = fits.open('carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
mymask = r.get_mask(hdu=hdu1)
intmask = mymask.astype(int)
fits.writeto('mask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits', intmask, hdu1.header, clobber=True)
#fits.writeto('northmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits', intmask, hdu1.header, clobber=True)
#fits.writeto('southmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits', intmask, hdu1.header, clobber=True)

###

region_name = "totfil_ga.reg"
r = pyregion.open(region_name)

hdu1 = fits.open('../emap_Orion_A_bw1.0.fits')[1]
mymask = r.get_mask(hdu=hdu1)
intmask = mymask.astype(int)
fits.writeto('mask_emap_Orion_A_bw1.0.fits', intmask, hdu1.header, clobber=True)
#fits.writeto('northmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits', intmask, hdu1.header, clobber=True)
#fits.writeto('southmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits', intmask, hdu1.header, clobber=True)
