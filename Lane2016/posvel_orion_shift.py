import numpy as np
import pyfits
import os
import sys
from astropy.io import fits
from astropy import wcs

name_cube = 'pv_mask_imfit_13co_pix_2_Tmb.fits'
name_out = 'pv_mask_imfit_13co_pix_2_Tmb_trans_shift.fits'
name_out_halfmax = 'pv_mask_imfit_13co_pix_2_Tmb_trans_shift_halfmax.fits'
name_cube = 'pv_mask_imfit_c18o_pix_2_Tmb.fits'
name_out = 'pv_mask_imfit_c18o_pix_2_Tmb_trans_shift.fits'
name_out_halfmax = 'pv_mask_imfit_c18o_pix_2_Tmb_trans_shift_halfmax.fits'

mmap=pyfits.open(name_cube)
cube=mmap[0].data
print 'cube shape',cube.shape

delt_dim1=mmap[0].header['cdelt1'] # position1
n_dim1=mmap[0].header['naxis1']
delt_dim2=mmap[0].header['cdelt2'] # velocity2
n_dim2=mmap[0].header['naxis2']
print 'cdelt1, cdelt2', delt_dim1, delt_dim2

v_step=np.abs(delt_dim2)
print 'v_step', v_step

new_vel_dim = int(6./v_step) # new number of pixels along velocity
vel_pos_map=np.zeros(shape=(n_dim1,new_vel_dim)) # new position 2, new velocity 1 
vel_pos_map_halfmax=np.zeros(shape=(n_dim1,new_vel_dim)) # new position 2, new velocity 1 

for index in range(n_dim1):
    peakind = np.argmax(cube[:,index])
    startind = int(peakind-3./v_step) # starting index in the old cube at peak velocity - 3 km/s
    vel_pos_map[index,:] = cube[startind:startind+new_vel_dim,index]
    greaterthanhalfmax = (vel_pos_map[index,:]>=0.5*cube[peakind,index])
    vel_pos_map_halfmax[index,greaterthanhalfmax] = 1

vel_pos_hdu=pyfits.PrimaryHDU(vel_pos_map)
vel_pos_hdu.header.update(CDELT1=mmap[0].header['CDELT2'])
vel_pos_hdu.header.update(CTYPE1='Velocity km/s')
vel_pos_hdu.header.update(NAXIS1=new_vel_dim)
vel_pos_hdu.header.update(CRVAL1=-3.)
vel_pos_hdu.header.update(CDELT2=mmap[0].header['CDELT1'])
vel_pos_hdu.header.update(CRVAL2=mmap[0].header['CRVAL1'])
vel_pos_hdu.header.update(NAXIS2=n_dim1)
vel_pos_hdu.header.update(CTYPE2='Offsets [arcmin]')
os.system('rm '+name_out)
vel_pos_hdu.writeto(name_out)

vel_pos_halfmax_hdu=pyfits.PrimaryHDU(vel_pos_map_halfmax)
vel_pos_halfmax_hdu.header.update(CDELT1=mmap[0].header['CDELT2'])
vel_pos_halfmax_hdu.header.update(CTYPE1='Velocity km/s')
vel_pos_halfmax_hdu.header.update(NAXIS1=new_vel_dim)
vel_pos_halfmax_hdu.header.update(CRVAL1=-3.)
vel_pos_halfmax_hdu.header.update(CDELT2=mmap[0].header['CDELT1'])
vel_pos_halfmax_hdu.header.update(CRVAL2=mmap[0].header['CRVAL1'])
vel_pos_halfmax_hdu.header.update(NAXIS2=n_dim1)
vel_pos_halfmax_hdu.header.update(CTYPE2='Offsets [arcmin]')
os.system('rm '+name_out_halfmax)
vel_pos_halfmax_hdu.writeto(name_out_halfmax)

