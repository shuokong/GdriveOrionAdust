import numpy as np
import math
import scipy.optimize as optimize
import matplotlib.pylab as plt

#import pyspeckit
import pyfits
import pylab
import aplpy
from numpy import ma
from matplotlib.pyplot import step, legend, xlim, ylim, show
import os
import radiomodule_orion as rmod
import posvel_orion as pv
#from matplotlib.backends.backend_pdf import PdfPages

def find_nearest(ralist,declist,pointra,pointdec,mindist):
    distances = ((ralist-pointra)**2+(declist-pointdec)**2)**0.5
    idx = distances.argmin()
    if distances[idx] < mindist:
        return idx
    else: 
        return -1

def main():

    name_cube = '/Users/shuokong/GoogleDrive/13co/products/mask_imfit_13co_pix_2_Tmb.fits'
    
    # one pv cut
    #raw_ra_list = np.array([83.52329236,83.57316692,83.62030445,83.6670613,83.71382045,83.75252857,83.79620141,83.82111418,83.84564682,83.86251383,83.86210348,83.85289673,83.84292146,83.82604083,83.80992492,83.80935538,83.80014058,83.79092437,83.77326557,83.7632777,83.74561138,83.73638688,83.73637222,83.74557319,83.75323795,83.77088618,83.79621755,83.82001665,83.84612392,83.8783804,83.91332917,83.95142695,83.98952916,84.0308125,84.07400417,84.11910176,84.1642031,84.21046505])
    #raw_dec_list = np.array([-4.876200218,-4.879461026,-4.897820323,-4.915411506,-4.93299871,-4.965857804,-4.99067018,-5.033833723,-5.076613716,-5.123974925,-5.173573767,-5.222459592,-5.271344718,-5.318700697,-5.36452826,-5.412955261,-5.461988054,-5.510106789,-5.555931851,-5.605577113,-5.6519109,-5.701555547,-5.752730272,-5.801156183,-5.84927762,-5.897401317,-5.940943148,-5.986011267,-6.029551074,-6.068507208,-6.103567778,-6.135131552,-6.166692778,-6.195094167,-6.219704167,-6.240522195,-6.261336664,-6.282571724]) # add one by hand to stretch out a bit more to the south, not included in pvpoints.reg
    # two pv cuts
    raw_ra_list = np.array([83.67620913,83.72689646,83.76744555,83.80672999,83.83207412,83.83206894,83.83206385,83.81430922,83.79655257,83.7800611,83.76229848,83.75467822,83.74578865,83.73816527,83.73815189,83.7279902,83.71887032,83.70120011,83.69197158,83.69195291,83.70115023,83.70881131,83.72645601,83.75178434,83.77558022,83.80168441,83.83393826,83.86888476,83.90698058,83.94508081,83.98636246,84.02955277,84.07464931,84.11974958])
    raw_dec_list = np.array([-4.917350655,-4.935047706,-4.962208694,-4.992523282,-5.036723747,-5.085971062,-5.133955069,-5.183198853,-5.227390751,-5.272844674,-5.320822798,-5.370065363,-5.418044445,-5.468548752,-5.517792235,-5.564507161,-5.61503437,-5.661366732,-5.711010624,-5.762185328,-5.810611952,-5.858733988,-5.906859078,-5.950402914,-5.995472923,-6.039014807,-6.077973505,-6.113036854,-6.144603656,-6.176167913,-6.204572583,-6.229186014,-6.25000762,-6.270825669]) 
    pdfname = 'pvKirkcores2.pdf'
    ttitle = 'PV cut 2'
    name_out = 'pv2_mask_imfit_13co_pix_2_Tmb.fits'
    #raw_ra_list = np.array([83.69015816,83.73672454,83.78361063,83.82954797,83.86883473,83.89418272,83.89418227,83.89418184,83.87643227,83.85868025,83.84219355,83.82443602,83.81682089,83.80793639,83.80031836,83.80031013,83.79015354,83.78103919,83.76337426,83.75415127,83.7541382,83.7633407,83.77100699,83.78865679,83.8139896,83.8377902,83.86389894,83.89615678,83.93110679,83.96920572,84.00730909,84.04859351,84.09178615,84.13688461,84.18198683])
    #raw_dec_list = np.array([-4.889574312,-4.909478974,-4.928116933,-4.954011799,-4.984322117,-5.028519878,-5.077767293,-5.125751391,-5.174997223,-5.219191156,-5.26464697,-5.312627124,-5.361870595,-5.409850716,-5.460355918,-5.509599449,-5.556315535,-5.606843787,-5.653178126,-5.702823061,-5.753997783,-5.802423402,-5.850544597,-5.898667737,-5.94220877,-5.98727614,-6.030815126,-6.069770245,-6.104829715,-6.136392291,-6.167952318,-6.196352408,-6.220961049,-6.241777657,-6.262590706]) 
    #pdfname = 'pvKirkcores1.pdf'
    #ttitle = 'PV cut 1'
    #name_out = 'pv1_mask_imfit_13co_pix_2_Tmb.fits'

    ra_list = np.concatenate([np.linspace(i,raw_ra_list[n+1],num=60,endpoint=False) for n,i in enumerate(raw_ra_list[:-1])])

    dec_list = np.concatenate([np.linspace(i,raw_dec_list[n+1],num=60,endpoint=False) for n,i in enumerate(raw_dec_list[:-1])])

    print 'len(ra_list) == len(dec_list)', len(ra_list) == len(dec_list), ra_list, dec_list

    vel_range = [5,14] #kms
    
    vel_rms = [-2,0.5] #kms
    
    mapa_maker = 'miriad'
    
    beam_param = [180.,180.,0.]
    
    center_coord = [83.806,-5.368]
    
    mindist = 1.5/60.
    
    corenames,xw,yw,coremasses,corevelocities12CO,coresnr12CO,coremom0_12CO,coremom1_12CO,coremom2_12CO,corevelocities13CO,coresnr13CO,coremom0_13CO,coremom1_13CO,coremom2_13CO,corevelocitiesC18O,coresnrC18O,coremom0_C18O,coremom1_C18O,coremom2_C18O,nh3velocities,nh3evelocities,nh3tkin,coregaus_x0_13CO,coregaus_ex0_13CO,coregaus_sigma_13CO,coregaus_esigma_13CO,coregaus_x0_C18O,coregaus_ex0_C18O,coregaus_sigma_C18O,coregaus_esigma_C18O = np.loadtxt('convol32_Kirkcores_velocities.txt',unpack=True)
    intcn = corenames.astype('int')
    
    #pvmap_orion,crval1,cdelt1 = pv.pos_vel(name_cube,name_out,ra_list,dec_list,vel_range,vel_rms,mapa_maker,beam_param,center_coord)
    hdu1 = pyfits.open(name_out)
    crval1 = hdu1[0].header['CRVAL1']
    cdelt1 = hdu1[0].header['CDELT1']
    hdu1.close()

    corepp = []
    coreind = []
    for cc,nn in enumerate(corenames):
        nearestpix = find_nearest(ra_list,dec_list,xw[cc],yw[cc],mindist)
        if nearestpix != -1:
            corepp.append(crval1+nearestpix*cdelt1) # assumes first element of ra_list is index 0 in array and ra_list follows the position direction
            coreind.append(cc)
    print 'min(corepp),max(corepp)', min(corepp),max(corepp)
    
    fig=plt.figure(figsize=(20,10))
    gc=aplpy.FITSFigure(name_out,dimensions=[0,1],figure=fig,hdu=0)
    gc.show_colorscale(aspect='auto')
    for ccc,nnn in enumerate(corepp):
        if nh3velocities[coreind[ccc]] < vel_range[0] or nh3velocities[coreind[ccc]] > vel_range[1] or np.isnan(nh3velocities[coreind[ccc]]): continue
        #gc.add_label(nnn,nh3velocities[coreind[ccc]],intcn[coreind[ccc]],color='blue')
        gc.show_markers(nnn,nh3velocities[coreind[ccc]],c='b')
    
    gc.ticks.set_xspacing(21.)
    gc.ticks.set_minor_frequency(7)
    gc.axis_labels.set_ytext('Velocity (km/s)')
    gc.ticks.show()
    plt.title(ttitle,size=20)
    gc.ticks.set_color('black')
    gc.ticks.set_length(10)
    gc.ticks.set_linewidth(2)
    gc.add_colorbar()
    gc.set_theme('publication')
    gc.colorbar.set_width(0.2)
    #gc.colorbar.set_axis_label_text('K')
    gc.colorbar.set_font(size=20)
    gc.set_tick_labels_font(size=20)
    gc.set_axis_labels_font(size=20)
    plt.savefig(pdfname,bbox_inches='tight')
    #plt.show()
    os.system('open '+pdfname)

if __name__=='__main__':
  main()

