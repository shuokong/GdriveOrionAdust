#!/usr/bin/python
# Filename: radiomodule.py


import numpy as np
import matplotlib.pylab as plt
import math
from matplotlib import rc


import pyfits
import pylab
#import aplpy
from numpy import ma
from matplotlib.pyplot import step, legend, xlim, ylim, show
import os
import sys
from astropy.io import fits
from astropy import wcs


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def coord_array(name_map,mapa_maker):
  map=fits.open(name_map)

  dim1=map[0].header['naxis1'] #ra
  dim2=map[0].header['naxis2'] #dec
  dim3=map[0].header['naxis3']
  naxis=map[0].header['naxis']
  ctype3=map[0].header['ctype3'] #'Vel-kms'
  print 'naxis', naxis


  w=wcs.WCS(map[0].header)


  dim2_arr=np.zeros(dim2)
  dim1_arr=np.zeros(dim1)
  dim3_arr=np.zeros(dim3)

  if naxis==4 :  # TODO: optimizar el for. Usar el len(dim), minimizar las ite
    for el in range(len(dim2_arr)):
      dim2_arr[el]=w.wcs_pix2world(0,el,0,0,0)[1]

    for el in range(len(dim1_arr)):
      dim1_arr[el]=w.wcs_pix2world(el,0,0,0,0)[0]

    for el in range(len(dim3_arr)):
      dim3_arr[el]=w.wcs_pix2world(0,0,el,0,0)[2]


  if mapa_maker=='miriad' or mapa_maker=='model' or mapa_maker=='complete':
    dim3_arr=dim3_arr*0.001





  map.close()
  return [dim1_arr,dim2_arr,dim3_arr]

def pos_to_spectrum(name_cube,coord_deg,mapa_maker,beam_param):

  map = pyfits.open(name_cube)
  if mapa_maker == 'stella' or mapa_maker=='model' or mapa_maker=='complete':
    cube = map[0].data
  else:
    cube = map[0].data[0]


  #print cube.shape
  [dim1_arr,dim2_arr,v] = coord_array(name_cube,mapa_maker)
  #print dim1_arr

  if len(beam_param)==0:
    ra=find_nearest(dim1_arr,coord_deg[0])
    dec=find_nearest(dim2_arr,coord_deg[1])

    spec=cube[:,dec,ra]
    return [spec,v]
  else:

    cosalpha=np.cos(beam_param[2]*np.pi/180)
    sinalpha=np.sin(beam_param[2]*np.pi/180)
    #print 'bpa', beam_param[2]
    a=beam_param[0]
    #print 'bmaj', a
    b=beam_param[1]
    #print 'bmin', b
    cont=0
    spec=np.zeros(len(v))

    ra_i=find_nearest(dim1_arr,coord_deg[0]+20./3600./np.cos(coord_deg[0]*np.pi/180.))
    ra_f=find_nearest(dim1_arr,coord_deg[0]-20./3600./np.cos(coord_deg[0]*np.pi/180.))
    dec_i=find_nearest(dim2_arr,coord_deg[1]-20./3600)
    dec_f=find_nearest(dim2_arr,coord_deg[1]+20./3600)


#    for ra in range(len(dim1_arr)):
    for ra in range(len(dim1_arr[ra_i:ra_f+1])):

#      for dec in range(len(dim2_arr)):
      for dec in range(len(dim2_arr[dec_i:dec_f+1])):

#        ra_offset=(dim1_arr[ra]-coord_deg[0])*np.cos(dim2_arr[dec]*np.pi/180.)*3600.
        ra_offset=(dim1_arr[ra+ra_i]-coord_deg[0])*np.cos(dim2_arr[dec+dec_i]*np.pi/180.)*3600.

#        dec_offset=(dim2_arr[dec]-coord_deg[1])*3600.
        dec_offset=(dim2_arr[dec+dec_i]-coord_deg[1])*3600.


        xp=ra_offset*cosalpha+dec_offset*sinalpha
        yp=-ra_offset*sinalpha+dec_offset*cosalpha


        radius=(xp/a)**2+(yp/b)**2

        if radius<=1:
          cont=cont+1
#          spec_t=cube[:,dec,ra]
          spec_t=cube[:,dec+dec_i,ra+ra_i]

          spec=spec+spec_t


          #gc2.show_markers(dim1_arr[ra],dim2_arr[dec],marker='+',s=200)

    #gc2.show_colorscale()

    #plt.show()
    #plt.plot(v,spec/cont)
    #plt.show()
    return [spec/cont,v]






