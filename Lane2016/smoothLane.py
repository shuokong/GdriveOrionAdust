import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

#cube = SpectralCube.read('beam_OrionA_850_auto_mos_clip.fits')
#cube.allow_huge_operations = True
#beam = radio_beam.Beam(major=60*u.arcsec, minor=60*u.arcsec, pa=0*u.deg)
#new_cube = cube.convolve_to(beam)

from astropy.io import fits
import astropy.wcs as wcs
from astropy import convolution

fh = fits.open('OrionA_850_auto_mos_clip.fits')
target_beam = radio_beam.Beam(major=60*u.arcsec, minor=60*u.arcsec, pa=0*u.deg)
orig_beam = radio_beam.Beam(14.6*u.arcsec) # I'm guessing....
kernel_beam = target_beam.deconvolve(orig_beam)
ww = wcs.WCS(fh[0].header)
pixscale = wcs.utils.proj_plane_pixel_scales(ww)[0] * u.deg
smooth_data = convolution.convolve_fft(fh[0].data.squeeze(), kernel_beam.as_kernel(pixscale), allow_huge=True)
fits.writeto('convol60_OrionA_850_auto_mos_clip.fits', smooth_data, fh[0].header, clobber=True)

