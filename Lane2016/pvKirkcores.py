from pvextractor import Path
from astropy import units as u
from astropy.coordinates import FK5
from pvextractor import extract_pv_slice
from spectral_cube import SpectralCube

makeslice = 0
if makeslice == 1:
    fk = FK5([83.52329236,83.57316692,83.62030445,83.6670613,83.71382045,83.75252857,83.79620141,83.82111418,83.84564682,83.86251383,83.86210348,83.85289673,83.84292146,83.82604083,83.80992492,83.80935538,83.80014058,83.79092437,83.77326557,83.7632777,83.74561138,83.73638688,83.73637222,83.74557319,83.75323795,83.77088618,83.79621755,83.82001665,83.84612392,83.8783804,83.91332917,83.95142695,83.98952916,84.0308125,84.07400417,84.11910176,84.1642031] * u.deg, [-4.876200218,-4.879461026,-4.897820323,-4.915411506,-4.93299871,-4.965857804,-4.99067018,-5.033833723,-5.076613716,-5.123974925,-5.173573767,-5.222459592,-5.271344718,-5.318700697,-5.36452826,-5.412955261,-5.461988054,-5.510106789,-5.555931851,-5.605577113,-5.6519109,-5.701555547,-5.752730272,-5.801156183,-5.84927762,-5.897401317,-5.940943148,-5.986011267,-6.029551074,-6.068507208,-6.103567778,-6.135131552,-6.166692778,-6.195094167,-6.219704167,-6.240522195,-6.261336664] * u.deg)
    path4 = Path(fk, width=6 * u.arcmin)
    slice3 = extract_pv_slice('/Users/shuokong/GoogleDrive/13co/products/regrid_12co_specsmooth_0p25_mask_imfit_13co_pix_2_Tmb.fits', path4)
    slice3.writeto('my_slice.fits')

import aplpy
import matplotlib.pyplot as plt
fig=plt.figure()
gc=aplpy.FITSFigure('my_slice.fits',dimensions=[0,1],figure=fig,hdu=0)
gc.show_colorscale(aspect='auto')

gc.ticks.set_xspacing(21.)
gc.ticks.set_minor_frequency(7)
gc.axis_labels.set_ytext('Velocity (km/s)')
gc.ticks.show()
gc.ticks.set_color('black')
gc.ticks.set_length(10)
gc.ticks.set_linewidth(2)
gc.add_colorbar()
gc.set_theme('publication')
gc.colorbar.set_width(0.2)
gc.colorbar.set_axis_label_text('K')
plt.savefig('pvKirkcores.pdf',bbox_inches='tight')

