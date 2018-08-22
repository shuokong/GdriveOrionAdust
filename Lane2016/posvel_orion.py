import numpy as np
import scipy.optimize as opt
import matplotlib; matplotlib.use("Agg")

import matplotlib.pylab as plt
import math

import pyfits
import pylab
from numpy import ma
from matplotlib.pyplot import step, legend, xlim, ylim, show
import os
import sys
from astropy.io import fits
from astropy import wcs
import radiomodule_orion as rmod


def pos_vel(name_cube,name_out,ra_list,dec_list,vel_range,vel_rms,mapa_maker,beam_param,center_coord):



  map=pyfits.open(name_cube)
  if mapa_maker=='miriad':
    cube=map[0].data[0]

  delt_dim1=map[0].header['cdelt1']
  c=np.abs(map[0].header['cdelt1'])
  delt_dim2=map[0].header['cdelt2']
  print 'cedelt1, cdelt2', delt_dim1*3600., delt_dim1*3600.

  [dim1_arr,dim2_arr,v]=rmod.coord_array(name_cube,mapa_maker)
  v_step=np.abs(v[0]-v[1])
  print 'v_step', v_step


  vel_i_ind=rmod.find_nearest(v,vel_range[0])
  vel_f_ind=rmod.find_nearest(v,vel_range[1])+1
  vel_i_rms=rmod.find_nearest(v,vel_rms[0])
  vel_f_rms=rmod.find_nearest(v,vel_rms[1])+1
  vel_plot=v[vel_i_ind:vel_f_ind]

  ra_listpos=ra_list
  dec_listpos=dec_list



  vel_pos_map=np.zeros(shape=(len(vel_plot),len(dec_listpos)))
  vel_pos_map.fill(np.nan)

  #rms=list()
  #maxmap=list()
  offsets_sum = 0. #np.zeros(len(ra_listpos))

  index_cent = rmod.find_nearest(dec_listpos,center_coord[1])
  print 'index of center coord', index_cent
  print 'center coordinate', ra_listpos[index_cent], dec_listpos[index_cent]

  for index in range(1,index_cent+1):
    offset = np.sqrt(((ra_listpos[index]-ra_listpos[index-1])*np.cos(dec_listpos[index-1]*np.pi/180.))**2+(dec_listpos[index]-dec_listpos[index-1])**2)
    offsets_sum = offsets_sum + offset
    if index == 1:
      delt = offset
      print 'offset delta', offset, delt

  print 'from upper edge to center the length is [arcmin] ', offsets_sum*60




  for index,elem in enumerate(ra_listpos):
      ra_pos=rmod.find_nearest(dim1_arr,elem)
      dec_pos=rmod.find_nearest(dim2_arr,dec_listpos[index])
      print 'index,vel_i_ind,vel_f_ind',index,vel_i_ind,vel_f_ind
      vel_pos_map[:,index]=np.array(rmod.pos_to_spectrum(name_cube,[elem,dec_listpos[index]],mapa_maker,beam_param))[0][vel_i_ind:vel_f_ind]
      #rms.append(np.std(np.array(rmod.pos_to_spectrum(name_cube,[elem,dec_listpos[index]],mapa_maker,beam_param))[0][vel_i_rms:vel_f_rms]))



  #rms=np.array(rms)

  #rms_val=np.mean(rms)




  vel_pos_hdu=pyfits.PrimaryHDU(vel_pos_map)
  vel_pos_hdu.header.update(NAXIS2=len(vel_plot))
  vel_pos_hdu.header.update(CDELT2=v_step)
  vel_pos_hdu.header.update(CRVAL2=vel_plot[0])#+8.0064)#+5.08688949)#+8.0064  #km/s #6.5-8.7) #
  vel_pos_hdu.header.update(CTYPE2='Velocity km/s')
  vel_pos_hdu.header.update(NAXIS1=len(dec_listpos))
  vel_pos_hdu.header.update(CDELT1=-delt*60)#!!!!!!!!!!!!
  #vel_pos_hdu.header.update('NAXIS2',len(vel_plot))
  #vel_pos_hdu.header.update('CDELT2',v_step)
  #vel_pos_hdu.header.update('CRVAL2',vel_plot[0])#+8.0064)#+5.08688949)#+8.0064  #km/s #6.5-8.7) #
  #vel_pos_hdu.header.update('CTYPE2','Velocity km/s')
  #vel_pos_hdu.header.update('NAXIS1',len(dec_listpos))
  #vel_pos_hdu.header.update('CDELT1',-delta*60)#!!!!!!!!!!!!
  #vel_pos_hdu.header.update('CDELT1',delt_dim2*np.pi/180*250*206265.)
  #vel_pos_hdu.header.update('CRVAL1',offsets_sum*60.)
  #vel_pos_hdu.header.update('CTYPE1','Offsets [arcmin]')

  vel_pos_hdu.header.update(CRVAL1=offsets_sum*60.)

  vel_pos_hdu.header.update(CTYPE1='Offsets [arcmin]')
  #if np.isnan(rms_val):
  #  rms_val=0.
  #  print 'rms is nan'
  #if np.isnan(maxmap_val):
  #  maxmap_val=0.
  #  print 'maxmap_val is nan'

  #vel_pos_hdu.header.update('rms',rms_val)
  #vel_pos_hdu.header.update('max',maxmap_val)

  os.system('rm '+name_out)
  vel_pos_hdu.writeto(name_out)



  # output file name, CRVAL1, CDELT1 assuming CRPIX1 = 1 
  return name_out,offsets_sum*60.,-delt*60
