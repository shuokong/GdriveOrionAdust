import math
import sys
import os
import numpy as np
from scipy.optimize import curve_fit
import telescope_sensitivity as ts
import pyfits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':30,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

def gaus(x,a,x0,sigma):
    return a/sigma*np.exp(-(x-x0)**2/(2*sigma**2))

def multigaus(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        x0 = params[i] # i-th line center
        a = params[i+1] # i-th line amplitude
        sigma = params[i+2] # i-th line sigma
        y = y + a * np.exp(-(x-x0)**2/(2.*sigma**2))
    return y

def cuberms(cube,velaxis=0): # given a ppv cube, where axes order vel, y, x, estimate cube rms in emission free channels
    stdplane = [(1,2),(0,2),(0,1)]
    planerms = np.nanstd(cube,axis=stdplane[velaxis]) # calculate standard deviation along y z axes, which will give std in each xplane 
    rms = np.nanmin(planerms[np.nonzero(planerms)]) # get rms for all velocity planes, the minimum rms is used as cube rms 
    print 'cube rms = ',rms
    return rms

def beampixel(bmaj,bmin,bpa,corecenter,cellsize,beamfraction=1.): # bmaj, bmin, cellsize in arcsec, corecenter = [pixelx, pixely], input bpa in degree
    pixellist = []
    rotation = float(bpa)/180.*np.pi
    cosa = np.cos(rotation)
    sina = np.sin(rotation)
    squareradius = int(bmaj/cellsize) # define a search domain first, make it twice bmaj
    xcenter,ycenter = corecenter
    semimajor = bmaj/2./cellsize*beamfraction**0.5
    semiminor = bmin/2./cellsize*beamfraction**0.5
    for x in range(xcenter-squareradius,xcenter+squareradius+1):
        for y in range(ycenter-squareradius,ycenter+squareradius+1):
            if ((x-xcenter)*cosa+(y-ycenter)*sina)**2/semimajor**2 + ((x-xcenter)*sina-(y-ycenter)*cosa)**2/semiminor**2 < 1.:
                pixellist.append([x,y])
    return pixellist

def mom0(deltav,intens):
    """calculate 0th moment of input vel weighted by intens, both should be 1d list or array"""
    return sum([deltav*intens[i] for i in range(len(intens))])

def mom1(rawvel,rawintens,thres):
    """calculate first moment of input vel weighted by intens, both should be 1d list or array"""
    if len(rawvel) != len(rawintens):
        print 'input have different size'
        return False
    usedata = (rawintens>thres)
    vel = rawvel[usedata]
    intens = rawintens[usedata]
    return sum([vel[i]*intens[i] for i in range(len(vel))])/sum(intens)

def mom2(rawvel,rawintens,thres):
    """calculate second moment of input vel weighted by intens, both should be 1d list or array"""
    if len(rawvel) != len(rawintens):
        print 'input have different size'
        return False
    usedata = (rawintens>thres)
    vel = rawvel[usedata]
    intens = rawintens[usedata]
    mmom1 = mom1(vel,intens,thres)
    return (sum([intens[i]*(vel[i]-mmom1)**2 for i in range(len(vel))])/sum(intens))**0.5

########################
#corenames, xw, yw = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores_degree.txt',usecols=(0,1,2),unpack=True)
corenames, xw, yw, coremasses, Tkin, eTkin, sTkin = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/GAScores.txt',usecols=(0,1,2,3,10,11,12),unpack=True)
#corenames, xw, yw, coremasses, Tkin, eTkin, sTkin = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/test.txt',usecols=(0,1,2,3,10,11,12),unpack=True)
worldcoord4 = np.stack((xw,yw,np.zeros_like(xw),np.zeros_like(xw)),axis=1)
#worldcoord3 = np.stack((xw,yw,np.zeros_like(xw)),axis=1) # Alyssa cubes

cellsize = 2. # voxel size in arcsec
#lines = ['/Users/shuokong/GoogleDrive/12co/products/mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/13co/products/mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/c18o/products/mask_imfit_c18o_pix_2_Tmb.fits']
lines = ['/Users/shuokong/GoogleDrive/12co/products/convol32_mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/13co/products/convol32_mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/c18o/products/convol32_mask_imfit_c18o_pix_2_Tmb.fits'] # convol32 to match GAS survey
#lines = ['/Users/shuokong/GoogleDrive/12co/products/convol32_mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/13co/products/convol32_specsmooth_0p25_convol_12co_mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/c18o/products/convol32_specsmooth_0p25_mask_imfit_c18o_pix_2_Tmb.fits'] # Alyssa cubes convolved to 32"
linebeams = []
linedata = []
linerms = []
linefreq = [115.27120180,110.20135430,109.78217340]
linejkfac = []
linexx = []
lineyy = []
line3rdaxis = []
print 'getting line cube metadata'
for j in range(len(lines)):
    print lines[j]
    hdulist = pyfits.open(lines[j])
    header = hdulist[0].header
    del header['HISTORY']
    w = wcs.WCS(header)
    scidata = hdulist[0].data[0,:,:,:] # remove stokes, usually index order: v, dec, ra
    pixcoord = w.all_world2pix(worldcoord4,1) # FITS standard uses 1
    #if j == 0:
    #    scidata = hdulist[0].data[0,:,:,:] # remove stokes, usually index order: v, dec, ra
    #    pixcoord = w.all_world2pix(worldcoord4,1) # FITS standard uses 1
    #else:
    #    scidata = hdulist[0].data[:,:,:] # usually index order: v, dec, ra
    #    pixcoord = w.all_world2pix(worldcoord3,1) # FITS standard uses 1
    n1,n2,n3 = scidata.shape
    linedata.append(scidata) 
    bmaj = header['BMAJ']*3600. # convert to arcsec
    bmin = header['BMIN']*3600.
    bpa = header['BPA']
    crval3 = header['CRVAL3']
    cdelt3 = header['CDELT3']
    crpix3 = header['CRPIX3']
    line3rdaxis.append([crval3, cdelt3, crpix3, n1])
    #linebeams.append([bmaj,bmin,bpa]) 
    #datarms = cuberms(scidata)
    #print 'datarms',datarms
    #linerms.append(datarms)
    #freq = linefreq[j]
    #jkfac = ts.jkelli(bmaj,bmin,freq)
    #linejkfac.append(jkfac) 
    xx = pixcoord[:,0]
    yy = pixcoord[:,1]
    linexx.append(xx)
    lineyy.append(yy)
    hdulist.close()
print 'finish getting line cube metadata'
#linebeams = [[10.010999813688, 8.091999962928, -12.8900003433], [7.620999868944001, 6.155000813304, 9.93999958038], [10.499000083668, 7.742001023136, -0.40000000596]]
#linerms = [1.1954995, 1.093392, 0.7671435]
linebeams = [[32.00000412762, 32.000000774868, -81.7325820923], [32.00000412762, 32.000000774868, 15.3498620987], [32.00000412762, 31.999997422116, 89.7622070312]] # convol32 to match GAS survey
linerms = [0.49438724, 0.19888465, 0.2835226] # convol32 to match GAS survey
#linebeams = [[32.00000412762, 32.000000774868, -81.7325820923], [32.000000774868, 31.999997422116, 77.4425811768], [32.00000412762, 31.999997422116, 89.7622070312]]
#linerms = [0.49438724, 0.11525758, 0.16264242]
print 'linebeams',linebeams
print 'linerms',linerms
#sys.exit()

worldcoord = np.stack((xw,yw),axis=1)
nh3file = '/Users/shuokong/GoogleDrive/AncillaryData/GBT/OrionA_Vlsr_DR1_rebase3_flag.fits'
print nh3file
hdulist = pyfits.open(nh3file)
header = hdulist[0].header
del header['HISTORY']
w = wcs.WCS(header)
scidata = hdulist[0].data[:,:] # in order: dec, ra
nh3data = scidata 
pixcoord = w.all_world2pix(worldcoord,1) # FITS standard uses 1
xx = pixcoord[:,0]
yy = pixcoord[:,1]
nh3xx = xx
nh3yy = yy
hdulist.close()

nh3efile = '/Users/shuokong/GoogleDrive/AncillaryData/GBT/OrionA_eVlsr_DR1_rebase3_flag.fits'
print nh3efile
hdulist = pyfits.open(nh3efile)
scidata = hdulist[0].data[:,:] # in order: dec, ra
enh3data = scidata 
hdulist.close()
print 'finish getting nh3 metadata'

singleGaussDict = {
                   '1':{'13CO':[11.15],'C18O':[11.26]},  
                   '2':{'13CO':[11.04],'C18O':[11.26]},  
                   '3':{'13CO':[9.28],'C18O':[9.28]},    
                   '4':{'13CO':[9.39],'C18O':[9.39]},    
                   '5':{'13CO':[9.5],'C18O':[9.39]},     
                   '6':{'13CO':[9.28],'C18O':[9.28]},    
                   '7':{'13CO':[9.5],'C18O':[9.28]},     
                   '8':{'13CO':[11.15],'C18O':[11.15]},  
                   '9':{'13CO':[11.04],'C18O':[11.26]},  
                  '10':{'13CO':[11.15],'C18O':[11.15]},  
                  '11':{'13CO':[10.49],'C18O':[10.6]},   
                  '12':{'13CO':[9.28],'C18O':[9.28]},    
                  '13':{'13CO':[11.26],'C18O':[11.37]},  
                  '14':{'13CO':[11.26],'C18O':[11.15]},  
                  '15':{'13CO':[10.71],'C18O':[10.93]},  
                  '16':{'13CO':[11.7],'C18O':[12.03]},   
                  '17':{'13CO':[10.93],'C18O':[10.93]},  
                  '18':{'13CO':[11.37],'C18O':[11.37]},  
                  '19':{'13CO':[8.18],'C18O':[7.96]},    
                  '20':{'13CO':[10.93],'C18O':[11.04]},  
                  '21':{'13CO':[8.29],'C18O':[8.07]},    
                  '22':{'13CO':[11.26],'C18O':[11.15]},  
                  '23':{'13CO':[11.04],'C18O':[11.15]},  
                  '24':{'13CO':[10.71],'C18O':[10.6]},   
                  '25':{'13CO':[11.15],'C18O':[11.37]},  
                  '26':{'13CO':[10.05],'C18O':[9.06]},   
                  '27':{'13CO':[10.71],'C18O':[10.6]},   
                  '28':{'13CO':[11.15],'C18O':[11.15]},  
                  '29':{'13CO':[11.59],'C18O':[11.81]},  
                  '30':{'13CO':[10.82],'C18O':[10.82]},  
                  '31':{'13CO':[7.74],'C18O':[7.74]},    
                  '32':{'13CO':[8.62],'C18O':[8.29]},    
                  '33':{'13CO':[8.84],'C18O':[8.62]},    
                  '34':{'13CO':[11.15],'C18O':[11.26]},  
                  '35':{'13CO':[11.15],'C18O':[11.04]},  
                  '36':{'13CO':[8.51],'C18O':[8.51]},    
                  '37':{'13CO':[10.82],'C18O':[11.15]},  
                  '38':{'13CO':[8.07],'C18O':[8.18]},    
                  '39':{'13CO':[11.15],'C18O':[11.26]},  
                  '40':{'13CO':[7.19],'C18O':[7.08]},    
                  '41':{'13CO':[10.82],'C18O':[11.37]},  
                  '42':{'13CO':[10.93],'C18O':[11.04]},  
                  '43':{'13CO':[10.49],'C18O':[10.6]},   
                  '44':{'13CO':[8.62],'C18O':[8.62]},    
                  '45':{'13CO':[8.29],'C18O':[8.07]},    
                  '46':{'13CO':[10.93],'C18O':[11.04]},  
                  '47':{'13CO':[8.51],'C18O':[8.51]},    
                  '48':{'13CO':[8.62],'C18O':[8.18]},    
                  '49':{'13CO':[10.71],'C18O':[10.6]},   
                  '50':{'13CO':[11.26],'C18O':[11.37]},  
                  '51':{'13CO':[11.26],'C18O':[11.37]},  
                  '52':{'13CO':[9.17],'C18O':[7.85]},    
                  '53':{'13CO':[10.49],'C18O':[10.49]},  
                  '54':{'13CO':[11.26],'C18O':[11.37]},  
                  '55':{'13CO':[11.26],'C18O':[11.26]},  
                  '56':{'13CO':[8.84],'C18O':[8.4]},     
                  '57':{'13CO':[11.04],'C18O':[10.93]},  
                  '58':{'13CO':[10.49],'C18O':[10.49]},  
                  '59':{'13CO':[8.62],'C18O':[8.07]},    
                  '60':{'13CO':[7.3],'C18O':[7.19]},     
                  '61':{'13CO':[11.37],'C18O':[11.48]},  
                  '62':{'13CO':[9.06],'C18O':[8.95]},    
                  '63':{'13CO':[9.61],'C18O':[9.5]},     
                  '64':{'13CO':[11.81],'C18O':[11.92]},  
                  '65':{'13CO':[11.7],'C18O':[11.59]},   
                  '66':{'13CO':[7.74],'C18O':[7.63]},    
                  '67':{'13CO':[7.96],'C18O':[7.96]},    
                  '68':{'13CO':[10.93],'C18O':[10.82]},  
                  '69':{'13CO':[10.93],'C18O':[10.93]},  
                  '70':{'13CO':[11.59],'C18O':[11.59]},  
                  '71':{'13CO':[10.82],'C18O':[10.93]},  
                  '72':{'13CO':[10.71],'C18O':[10.82]},  
                  '73':{'13CO':[6.86],'C18O':[6.86]},    
                  '74':{'13CO':[8.4],'C18O':[8.51]},     
                  '75':{'13CO':[11.7],'C18O':[11.81]},   
                  '76':{'13CO':[9.5],'C18O':[9.17]},     
                  '77':{'13CO':[7.85],'C18O':[7.74]},    
                  '78':{'13CO':[10.93],'C18O':[10.82]},  
                  '79':{'13CO':[10.93],'C18O':[10.93]},  
                  '80':{'13CO':[8.95],'C18O':[8.73]},    
                  '81':{'13CO':[10.6],'C18O':[10.6]},    
                  '82':{'13CO':[6.75],'C18O':[6.75]},    
                  '83':{'13CO':[7.3],'C18O':[7.08]},     
                  '84':{'13CO':[11.15],'C18O':[11.37]},  
                  '85':{'13CO':[10.49],'C18O':[10.71]},  
                  '86':{'13CO':[10.49],'C18O':[10.6]},   
                  '87':{'13CO':[11.59],'C18O':[11.7]},   
                  '88':{'13CO':[11.15],'C18O':[11.26]},  
                  '89':{'13CO':[8.29],'C18O':[8.29]},    
                  '90':{'13CO':[8.84],'C18O':[8.73]},    
                  '91':{'13CO':[11.26],'C18O':[11.37]},  
                  '92':{'13CO':[10.6],'C18O':[10.6]},    
                  '93':{'13CO':[7.85],'C18O':[7.96]},    
                  '94':{'13CO':[10.49],'C18O':[10.49]},  
                  '95':{'13CO':[11.26],'C18O':[11.37]},  
                  '96':{'13CO':[7.96],'C18O':[7.96]},    
                  '97':{'13CO':[10.82],'C18O':[11.04]},  
                  '98':{'13CO':[11.15],'C18O':[11.04]},  
                  '99':{'13CO':[11.04],'C18O':[11.15]},  
                 '100':{'13CO':[10.82],'C18O':[10.93]},  
                 '101':{'13CO':[11.59],'C18O':[11.59]},  
                 '102':{'13CO':[7.19],'C18O':[7.19]},    
                 '103':{'13CO':[8.62],'C18O':[8.51]},    
                 '104':{'13CO':[10.82],'C18O':[10.93]},  
                 '105':{'13CO':[8.62],'C18O':[8.18]},    
                 '106':{'13CO':[8.73],'C18O':[8.73]},    
                 '107':{'13CO':[8.51],'C18O':[8.4]},     
                 '108':{'13CO':[11.48],'C18O':[11.26]},  
                 '109':{'13CO':[12.36],'C18O':[12.47]},  
                 '110':{'13CO':[7.96],'C18O':[7.74]},    
                 '111':{'13CO':[7.52],'C18O':[7.41]},    
                 '112':{'13CO':[8.62],'C18O':[8.73]},    
                 '113':{'13CO':[10.93],'C18O':[10.93]},  
                 '114':{'13CO':[10.93],'C18O':[10.93]},  
                 '115':{'13CO':[10.49],'C18O':[11.26]},  
                 '116':{'13CO':[7.74],'C18O':[7.74]},    
                 '117':{'13CO':[10.27],'C18O':[10.05]},  
                 '118':{'13CO':[10.93],'C18O':[10.71]},  
                 '119':{'13CO':[11.7],'C18O':[12.03]},   
                 '120':{'13CO':[9.28],'C18O':[9.39]},    
                 '121':{'13CO':[9.28],'C18O':[9.39]},    
                 '122':{'13CO':[9.94],'C18O':[9.61]},    
                 '123':{'13CO':[7.63],'C18O':[7.74]},    
                 '124':{'13CO':[10.27],'C18O':[10.27]},  
                 '125':{'13CO':[8.29],'C18O':[8.18]},    
                 '126':{'13CO':[8.51],'C18O':[8.29]},    
                 '127':{'13CO':[12.14],'C18O':[12.36]},  
                 '128':{'13CO':[11.26],'C18O':[11.37]},  
                 '129':{'13CO':[8.18],'C18O':[8.84]},    
                 '130':{'13CO':[8.73],'C18O':[8.07]},    
                 '131':{'13CO':[7.96],'C18O':[7.85]},    
                 '132':{'13CO':[7.19],'C18O':[6.97]},    
                 '133':{'13CO':[9.39],'C18O':[9.61]},    
                 '134':{'13CO':[11.04],'C18O':[11.04]},  
                 '135':{'13CO':[11.59],'C18O':[11.7]},   
                 '136':{'13CO':[8.62],'C18O':[8.29]},    
                 '137':{'13CO':[11.04],'C18O':[11.04]},  
                 '138':{'13CO':[10.93],'C18O':[10.49]},  
                 '139':{'13CO':[9.39],'C18O':[9.28]},    
                 '140':{'13CO':[11.92],'C18O':[10.6]},   
                 '141':{'13CO':[7.85],'C18O':[7.63]},    
                 '142':{'13CO':[7.96],'C18O':[8.07]},    
                 '143':{'13CO':[10.82],'C18O':[10.71]},  
                 '144':{'13CO':[9.83],'C18O':[9.61]},    
                 '145':{'13CO':[11.26],'C18O':[11.26]},  
                 '146':{'13CO':[10.27],'C18O':[10.38]},  
                 '147':{'13CO':[10.38],'C18O':[10.05]},  
                 '148':{'13CO':[8.18],'C18O':[7.41]},    
                 '149':{'13CO':[9.17],'C18O':[8.95]},    
                 '150':{'13CO':[11.37],'C18O':[11.59]},  
                 '151':{'13CO':[8.62],'C18O':[8.51]},    
                 '152':{'13CO':[9.17],'C18O':[8.95]},    
                 '153':{'13CO':[9.06],'C18O':[9.06]},    
                 '154':{'13CO':[6.64],'C18O':[6.53]},    
                 '155':{'13CO':[11.04],'C18O':[11.15]},  
                 '156':{'13CO':[7.41],'C18O':[8.84]},    
                 '157':{'13CO':[8.84],'C18O':[7.19]},    
                 '158':{'13CO':[11.48],'C18O':[11.59]},  
                 '159':{'13CO':[10.71],'C18O':[10.93]},  
                 '160':{'13CO':[10.6],'C18O':[10.6]},    
                 '161':{'13CO':[10.93],'C18O':[10.93]},  
                 '162':{'13CO':[7.96],'C18O':[8.95]},    
                 '163':{'13CO':[9.72],'C18O':[9.5]},     
                 '164':{'13CO':[7.08],'C18O':[6.97]},    
                 '165':{'13CO':[8.07],'C18O':[8.62]},    
                 '166':{'13CO':[10.82],'C18O':[10.71]},  
                 '167':{'13CO':[8.84],'C18O':[8.73]},    
                 '168':{'13CO':[8.18],'C18O':[7.41]},    
                 '169':{'13CO':[11.26],'C18O':[11.26]},  
                 '170':{'13CO':[8.29],'C18O':[8.07]},    
                 '171':{'13CO':[7.85],'C18O':[9.06]},    
                 '172':{'13CO':[10.71],'C18O':[10.82]},  
                 '173':{'13CO':[11.7],'C18O':[11.7]},    
                 '174':{'13CO':[10.71],'C18O':[10.49]},  
                 '175':{'13CO':[11.04],'C18O':[11.04]},  
                 '176':{'13CO':[7.52],'C18O':[7.41]},    
                 '177':{'13CO':[8.4],'C18O':[8.62]},     
                 '178':{'13CO':[6.97],'C18O':[7.08]},    
                 '179':{'13CO':[8.18],'C18O':[6.86]},    
                 '180':{'13CO':[11.04],'C18O':[11.04]},  
                 '181':{'13CO':[8.62],'C18O':[8.29]},    
                 '182':{'13CO':[9.72],'C18O':[9.5]},     
                 '183':{'13CO':[11.48],'C18O':[11.59]},  
                 '184':{'13CO':[9.06],'C18O':[9.17]},    
                 '185':{'13CO':[10.38],'C18O':[10.27]},  
                 '186':{'13CO':[7.52],'C18O':[7.41]},    
                 '187':{'13CO':[9.28],'C18O':[7.85]},    
                 '188':{'13CO':[7.96],'C18O':[7.96]},    
                 '189':{'13CO':[8.95],'C18O':[7.63]},    
                 '190':{'13CO':[7.3],'C18O':[7.3]},      
                 '191':{'13CO':[11.7],'C18O':[11.7]},    
                 '192':{'13CO':[11.48],'C18O':[11.48]},  
                 '193':{'13CO':[11.04],'C18O':[11.04]},  
                 '194':{'13CO':[11.15],'C18O':[11.04]},  
                 '195':{'13CO':[7.85],'C18O':[7.74]},    
                 '196':{'13CO':[8.84],'C18O':[8.73]},    
                 '197':{'13CO':[8.18],'C18O':[8.07]},    
                 '198':{'13CO':[8.51],'C18O':[8.51]},    
                 '199':{'13CO':[7.85],'C18O':[7.85]},    
                 '200':{'13CO':[7.63],'C18O':[7.19]},    
                 '201':{'13CO':[8.18],'C18O':[7.96]},    
                 '202':{'13CO':[11.15],'C18O':[11.04]},  
                 '203':{'13CO':[7.19],'C18O':[7.08]},    
                 '204':{'13CO':[9.17],'C18O':[9.17]},    
                 '205':{'13CO':[7.85],'C18O':[7.74]},    
                 '206':{'13CO':[9.5],'C18O':[7.96]},     
                 '207':{'13CO':[10.93],'C18O':[10.93]},  
                 '208':{'13CO':[8.51],'C18O':[7.96]},    
                 '209':{'13CO':[8.4],'C18O':[7.96]},     
                 '210':{'13CO':[9.06],'C18O':[8.95]},    
                 '211':{'13CO':[10.6],'C18O':[10.38]},   
                 '212':{'13CO':[8.73],'C18O':[8.29]},    
                 '213':{'13CO':[8.95],'C18O':[7.74]},    
                 '214':{'13CO':[10.49],'C18O':[10.38]},  
                 '215':{'13CO':[9.28],'C18O':[7.85]},    
                 '216':{'13CO':[10.27],'C18O':[10.49]},  
                 '217':{'13CO':[10.71],'C18O':[11.04]},  
                 '218':{'13CO':[8.62],'C18O':[7.08]},    
                 '219':{'13CO':[12.14],'C18O':[12.14]},  
                 '220':{'13CO':[11.15],'C18O':[11.15]},  
                 '221':{'13CO':[7.85],'C18O':[7.41]},    
                 '222':{'13CO':[11.26],'C18O':[11.26]},  
                 '223':{'13CO':[8.51],'C18O':[7.19]},    
                 '224':{'13CO':[7.96],'C18O':[7.96]},    
                 '225':{'13CO':[7.19],'C18O':[7.19]},    
                 '226':{'13CO':[8.29],'C18O':[8.51]},    
                 '227':{'13CO':[8.4],'C18O':[8.4]},      
                 '228':{'13CO':[8.62],'C18O':[8.07]},    
                 '229':{'13CO':[8.73],'C18O':[8.62]},    
                 '230':{'13CO':[11.7],'C18O':[11.81]},   
                 '231':{'13CO':[9.39],'C18O':[9.39]},    
                 '232':{'13CO':[9.5],'C18O':[8.18]},     
                 '233':{'13CO':[6.53],'C18O':[6.42]},    
                 '234':{'13CO':[8.62],'C18O':[7.3]},     
                 '235':{'13CO':[11.59],'C18O':[11.48]},  
                 '236':{'13CO':[9.39],'C18O':[9.5]},     
                 '237':{'13CO':[9.06],'C18O':[8.18]},    
                  }

multiGaussDict = {
                  '41':{'13CO':[8,10],'C18O':[8,10]},
                  '43':{'13CO':[8,9.5],'C18O':[8,9.5]},
                  '54':{'13CO':[9,10.5],'C18O':[9,10.5]},
                  '59':{'13CO':[8,9,10],'C18O':[8]},
                  '65':{'13CO':[9,11],'C18O':[9,11]},
                  '70':{'13CO':[8],'C18O':[8]}, # outflow? # 13CO bad fit
                  '73':{'13CO':[8],'C18O':[8]}, # outflow?
                  '77':{'13CO':[11,13],'C18O':[11,13]},
                  '81':{'13CO':[7,10.5],'C18O':[7]},
                  '83':{'13CO':[11,13],'C18O':[11]}, # 13CO bad fit
                  '97':{'13CO':[7,8.5],'C18O':[7,8.5]},
                  '99':{'13CO':[6,8.5],'C18O':[8.5]},
                  '107':{'13CO':[8,9],'C18O':[8,9]},
                  '110':{'13CO':[11.2,13.2],'C18O':[11,13]}, # 13CO bad fit
                  '114':{'13CO':[8.5,11],'C18O':[8.5,11]},
                  '117':{'13CO':[9.1,10.5,13],'C18O':[10.5]}, # 13CO bad fit
                  '121':{'13CO':[7,8.5],'C18O':[8.5]},
                  '122':{'13CO':[7,10.5],'C18O':[7]},
                  '127':{'13CO':[9,10.5],'C18O':[9,10.5]},
                  '133':{'13CO':[6.5,8],'C18O':[6.3,7.9]},
                  '144':{'13CO':[9,11],'C18O':[9,11]},
                  '153':{'13CO':[7,8],'C18O':[7,8]},
                  '154':{'13CO':[8],'C18O':[8]}, # outflow? # C18O bad fit # 13CO bad fit
                  '165':{'13CO':[7.5,9.5],'C18O':[7.5,9.5]},
                  '169':{'13CO':[9,11],'C18O':[9,11]},
                  '172':{'13CO':[6.5,8.8,10.5],'C18O':[7]}, # 13CO bad fit
                  '174':{'13CO':[7.2,10.7],'C18O':[7,8]},
                  '180':{'13CO':[8.5,10.5],'C18O':[8.5,10.5]},
                  '182':{'13CO':[11,13],'C18O':[11,13]},
                  '194':{'13CO':[8,10],'C18O':[8]},
                  '198':{'13CO':[10.3,11.3,13],'C18O':[11]},
                  '199':{'13CO':[8,9,10],'C18O':[8]}, # 13CO bad fit
                  '203':{'13CO':[8.5,9.5,11],'C18O':[9,11]},
                  '207':{'13CO':[11,13],'C18O':[11,13]},
                  '208':{'13CO':[11,13],'C18O':[11,13]},
                  '220':{'13CO':[7,8,9],'C18O':[7,8,9]}, # 13CO bad fit
                  '223':{'13CO':[10.8,12.8],'C18O':[11]}, # 13CO bad fit
                  '226':{'13CO':[7,9],'C18O':[7,9]},
                  '235':{'13CO':[8,10],'C18O':[8,10]},
                  '240':{'13CO':[11.0,13.0],'C18O':[11]},
                  '241':{'13CO':[11,13],'C18O':[11,13]},
                  '254':{'13CO':[7,9.5],'C18O':[7,9.5]}, # C18O bad fit # 13CO bad fit
                  '264':{'13CO':[8,11],'C18O':[8]},
                  '265':{'13CO':[6.5,8],'C18O':[8]},
                  '274':{'13CO':[8,9],'C18O':[9]},
                  '282':{'13CO':[6,8,9],'C18O':[6,8,9]},
                  '287':{'13CO':[7,10],'C18O':[7]},
                  '294':{'13CO':[6,8],'C18O':[8]},
                  '310':{'13CO':[10.5,12],'C18O':[10.5]},
                  '312':{'13CO':[5,8,11],'C18O':[8]},
                  '317':{'13CO':[8,10],'C18O':[8,10]},
                  '357':{'13CO':[6,8.5],'C18O':[8.5]},
                  '361':{'13CO':[8,9],'C18O':[8,9]},
                  '362':{'13CO':[7,9],'C18O':[7,9]},
                  '363':{'13CO':[7,10],'C18O':[7]},
                  '371':{'13CO':[7,9],'C18O':[7.3,8.8]},
                  '387':{'13CO':[7,9],'C18O':[7,9]},
                  '401':{'13CO':[10.5,11.5],'C18O':[11]},
                  '405':{'13CO':[7,8,9],'C18O':[7,8,9]},
                  '407':{'13CO':[9.5,11],'C18O':[9.5,11]},
                  '413':{'13CO':[7,10],'C18O':[7,10]},
                  '420':{'13CO':[7,8.5],'C18O':[7,8.5]},
                  '431':{'13CO':[7,9],'C18O':[7,9]},
                  '441':{'13CO':[8,9],'C18O':[7,8,9]},
                  '464':{'13CO':[9.5,11,13],'C18O':[11,13]},
                  '469':{'13CO':[8,10],'C18O':[7.5]},
                  '470':{'13CO':[7,8.5,9.5],'C18O':[7,8.5,9.5]},
                  '478':{'13CO':[7,8.5],'C18O':[7,8.5]}, # C18O bad fit
                  '493':{'13CO':[6,8],'C18O':[6,8]},
                  '498':{'13CO':[7,10],'C18O':[10]},
                  '508':{'13CO':[6.5,8,9],'C18O':[6.5,8,9]},
                  '524':{'13CO':[7.5,9,11],'C18O':[7.5,9]},
                  '535':{'13CO':[7,8,9],'C18O':[7,8,9]},
                  '538':{'13CO':[7.5,9],'C18O':[7.5,9]},
                  '539':{'13CO':[7.5,12],'C18O':[7.5]},
                  '545':{'13CO':[11.0,12.8],'C18O':[11]},
                  '570':{'13CO':[6,8],'C18O':[8]},
                  '572':{'13CO':[8,9],'C18O':[8,9]},
                  '584':{'13CO':[5.8,7.5],'C18O':[7]}, # 13CO bad fit
                  '586':{'13CO':[8.2,10.1],'C18O':[8.0,8.9]}, # C18O bad fit # 13CO bad fit
                  '592':{'13CO':[5.5,7,8.5],'C18O':[5.5,7,8.5]}, # 13CO bad fit
                  '594':{'13CO':[8.5,9.5],'C18O':[8.5,9.5]},
                  '595':{'13CO':[8,9],'C18O':[8,9]},
                  '601':{'13CO':[8,9],'C18O':[8,9]},
                  '610':{'13CO':[6.5,8],'C18O':[8]},
                  '611':{'13CO':[8,10],'C18O':[7.8,9.5]},
                  '623':{'13CO':[8,9],'C18O':[8,9]},
                  '632':{'13CO':[7,8,9],'C18O':[7,8,9]},
                  '637':{'13CO':[8.5,10.5],'C18O':[8.5,10.5]},
                  '648':{'13CO':[9,11],'C18O':[11]},
                  '686':{'13CO':[7.0,8.8,9.9],'C18O':[7.1,8.8]}, # C18O bad fit
                  '691':{'13CO':[10,12],'C18O':[10,12]},
                  '701':{'13CO':[7.5,8.5,10],'C18O':[8]},
                  '707':{'13CO':[7,8.5,10],'C18O':[7.1,8.6]},
                  '721':{'13CO':[7.9,9.5],'C18O':[7.9,9.5]}, # 13CO bad fit
                  '734':{'13CO':[7,9],'C18O':[7,8.5]},
                  '750':{'13CO':[7,8.5,11],'C18O':[7,8.5]},
                  '835':{'13CO':[8,9.5],'C18O':[8,9.5]},
                  '869':{'13CO':[6.5,8,10],'C18O':[6.5]},
                  '871':{'13CO':[7,9],'C18O':[7,9]},
                  '879':{'13CO':[8,9.5],'C18O':[8,9.5]},
                  '907':{'13CO':[8,9,10],'C18O':[8]},
                  }

#single_components = [str(int(iii)) for iii in corenames if str(int(iii)) not in multiGaussDict.keys()]
#print 'single_components',single_components,'len(single_components)',len(single_components),'len(multiGaussDict.keys())',len(multiGaussDict.keys()),set(single_components).intersection(set(multiGaussDict.keys()))
#outlier = [iii for iii in multiGaussDict.keys() if float(iii) not in corenames]
#print 'outlier',outlier
#sys.exit()
#single_components = [5, 9, 10, 12, 16, 17, 18, 19, 21, 23, 30, 32, 33, 34, 35, 36, 39, 40, 42, 44, 47, 49, 50, 55, 56, 57, 58, 63, 68, 69, 72, 82, 84, 91, 93, 95, 100, 101, 105, 108, 109, 116, 125, 128, 129, 130, 134, 137, 139, 142, 147, 161, 164, 166, 167, 170, 177, 179, 181, 189, 191, 192, 193, 197, 211, 213, 221, 225, 228, 231, 232, 234, 239, 245, 246, 247, 249, 250, 256, 257, 259, 260, 266, 269, 283, 291, 292, 293, 300, 302, 306, 319, 323, 327, 330, 331, 342, 353, 355, 367, 390, 391, 402, 422, 432, 434, 440, 461, 462, 463, 471, 486, 498, 500, 518, 537, 540, 543, 546, 551, 554, 571, 588, 603, 612, 616, 622, 630, 699, 702, 755, 768, 784, 790, 816, 873]

linenames = [r'$\rm ^{12}CO(1$-$0)$',r'$\rm ^{13}CO(1$-$0)$',r'$\rm C^{18}O(1$-$0)$']
molnames = ['12CO','13CO','C18O']
xpanels = 1
ypanels = len(lines)
xpanelwidth = 12
ypanelwidth = 5
vlow = 0
vhigh = 16
cols = len(corenames)
corevelocities = [[],[],[]]
coremom0s = [[],[],[]]
coremom1s = [[],[],[]]
coremom2s = [[],[],[]]
nh3velocities = [[]]
nh3evelocities = [[]]
nh3tkin = [[]]
coresnr = [[],[],[]]
coregaus_x0 = [[],[],[]]
coregaus_ex0 = [[],[],[]]
coregaus_sigma = [[],[],[]]
coregaus_esigma = [[],[],[]]
for cc in range(cols):
    #if str(int(corenames[cc])) not in multiGaussDict.keys(): continue
    fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels/1.1))
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    #pdfname='corespectra/Kirk/single_convol32_averspec_Kirk_core'+str(cc+1)+'.pdf'
    pdfname='corespectra/Kirk/convol32_averspec_Kirk_core'+str(cc+1)+'.pdf'
    #pdfname='corespectra/Kirk/convol32_averspec_all0p25channel_core'+str(cc+1)+'.pdf'
    datafiles = {}
    #print 'x,y',x,y
    ccnh3xx = int(nh3xx[cc])
    ccnh3yy = int(nh3yy[cc])
    try:
        ccvlsr = [nh3data[ccnh3yy,ccnh3xx]]
        ccevlsr = [enh3data[ccnh3yy,ccnh3xx]]
    except:
        ccvlsr = [-100]
        ccevlsr = [-100]
    print 'ccnh3xx,ccnh3yy',ccnh3xx,ccnh3yy,'ccvlsr',ccvlsr,'ccevlsr',ccevlsr
    nh3velocities[0].append(ccvlsr[0])
    nh3evelocities[0].append(ccevlsr[0])
    nh3tkin[0].append(Tkin[cc])
    ccvlsrstyle = ['dashed']
    for i in range(0,xpanels):
        for j in range(0,ypanels):
            print lines[j]
            panel = i+j*xpanels+1
            print 'panel',panel 
            bmaj,bmin,bpa = linebeams[panel-1]
            xx = linexx[panel-1]
            yy = lineyy[panel-1]
            corecenter = [int(xx[cc]),int(yy[cc])]
            pixlist = beampixel(bmaj,bmin,bpa+90.,corecenter,cellsize,beamfraction=1.)
            corei = np.array(pixlist) 
            data = linedata[panel-1] 
            try:
                temp = data[:,corei[:,1],corei[:,0]] 
            except:
                corevelocities[panel-1].append(-100)
                coresnr[panel-1].append(0)
                coremom0s[panel-1].append(-100)
                coremom1s[panel-1].append(-100)
                coremom2s[panel-1].append(-100)
                continue
            rawspectrum = np.nanmean(temp,axis=(1))
            crval3, cdelt3, crpix3, n1 = line3rdaxis[panel-1]
            rawvelocity = np.array([(crval3+cdelt3*((ii+1)-crpix3))/1.e3 for ii in range(n1)]) # should not use i, it will confuse with the for i above
            subvel = (rawvelocity > vlow) & (rawvelocity < vhigh) 
            rawintens = rawspectrum[subvel]
            velocity = rawvelocity[subvel]
            #coresavetxtarr = np.stack((velocity,rawintens),axis=1)
            #np.savetxt('corespectra/datapoints/'+molnames[j]+'/convol32_Kirkcore'+str(int(corenames[cc]))+'_'+molnames[j]+'.txt',coresavetxtarr,fmt='%7.2f %10.2e')
            coremom0 = mom0(cdelt3/1.e3,rawintens)
            coremom0s[panel-1].append(coremom0)
            velocityguess = []
            if j > 0:
                if str(int(corenames[cc])) in multiGaussDict.keys():
                    guess = []
                    for ii in multiGaussDict[str(int(corenames[cc]))][molnames[j]]:
                        guess += [ii,0.2,0.2]
                        velocityguess += [ii]
                    popt, pcov = curve_fit(multigaus, velocity, rawintens, p0=guess, maxfev=100000)
                    perr = np.sqrt(np.diag(pcov))
                    component_number = len(multiGaussDict[str(int(corenames[cc]))][molnames[j]])
                    component_velocities = [abs(popt[ii*3]-ccvlsr[0]) for ii in range(component_number)]
                    closest_index = np.argmin(component_velocities)
                    closest_velocity = popt[closest_index*3]
                    closest_evelocity = perr[closest_index*3]
                    closest_sigma = popt[closest_index*3+2]
                    closest_esigma = perr[closest_index*3+2]
                    coregaus_x0[panel-1].append(closest_velocity)
                    coregaus_ex0[panel-1].append(closest_evelocity)
                    coregaus_sigma[panel-1].append(closest_sigma)
                    coregaus_esigma[panel-1].append(closest_esigma)
                else:
                    guess = []
                    for ii in singleGaussDict[str(cc+1)][molnames[j]]:
                        guess += [ii,0.5,0.5]
                        velocityguess += [ii]
                    popt, pcov = curve_fit(multigaus, velocity, rawintens, p0=guess)
                    perr = np.sqrt(np.diag(pcov))
                    closest_velocity = popt[0]
                    closest_evelocity = perr[0]
                    closest_sigma = popt[2]
                    closest_esigma = perr[2]
                    coregaus_x0[panel-1].append(closest_velocity)
                    coregaus_ex0[panel-1].append(closest_evelocity)
                    coregaus_sigma[panel-1].append(closest_sigma)
                    coregaus_esigma[panel-1].append(closest_esigma)
            datarms = linerms[panel-1]
            threshold = datarms*5.
            try:
                coremom1 = mom1(velocity,rawintens,threshold)
                coremom1s[panel-1].append(coremom1)
            except:
                coremom1s[panel-1].append(-100)
            try:
                coremom2 = mom2(velocity,rawintens,threshold)
                coremom2s[panel-1].append(coremom2)
            except:
                coremom2s[panel-1].append(-100)
            print 'coremom2',lines[j],coremom2
            try:
                peakind = np.nanargmax(rawintens)
            except:
                corevelocities[panel-1].append(-100)
                coresnr[panel-1].append(0)
                continue
            corevelocities[panel-1].append(velocity[peakind])
            snr = rawintens[peakind]/datarms
            coresnr[panel-1].append(snr)
            ymin = 0
            ymax = 0
            if np.nanmin(rawintens) < ymin: ymin = np.nanmin(rawintens)
            if np.nanmax(rawintens) > ymax: ymax = np.nanmax(rawintens)
            datafiles['panel'+str(panel)] = {'title':linenames[j],'lines':{'1':{'x':velocity,'y':rawintens,'peakvelocity':velocity[peakind],'peaksnr':snr,'legends':'data','linestyles':'k-','drawsty':'steps-mid'},},'xlim':[vlow,vhigh],'ylim':[ymin-(ymax-ymin)/10.,ymax+(ymax-ymin)/10.],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm LSR}~\rm (km~s^{-1})$','ylabel':r'$T_{\rm mb}~\rm (K)$','text':'','vertlines':[ccvlsr[0]],'vertlinestyles':[ccvlsrstyle[0]],'vertlinecolors':['b'],'vertlinewidths':[4],'vertlinelengths':[0.5]}
            #datafiles['panel'+str(panel)] = {'title':linenames[j],'lines':{'1':{'x':velocity,'y':rawintens,'peakvelocity':velocity[peakind],'peaksnr':snr,'legends':'data','linestyles':'k-','drawsty':'steps-mid'},},'xlim':[vlow,vhigh],'ylim':[ymin-(ymax-ymin)/10.,ymax+(ymax-ymin)/10.],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm LSR}~\rm (km~s^{-1})$','ylabel':r'$T_{\rm mb}~\rm (K)$','text':'','vertlines':[ccvlsr[0],coremom1]+velocityguess,'vertlinestyles':[ccvlsrstyle[0],'dotted']+['solid' for ii in velocityguess],'vertlinecolors':['b','b']+['g' for ii in velocityguess],'vertlinewidths':[4,2]+[1 for ii in velocityguess],'vertlinelengths':[0.5,0.6]+[0.3 for ii in velocityguess]}
            if j > 0:
                datafiles['panel'+str(panel)]['lines']['2']={'x':velocity,'y':multigaus(velocity,*popt),'legends':'fit','linestyles':'g-','drawsty':'default'}
                datafiles['panel'+str(panel)]['vertlines'].append(closest_velocity)
                datafiles['panel'+str(panel)]['vertlinestyles'].append('dashed')
                datafiles['panel'+str(panel)]['vertlinecolors'].append('g')
                datafiles['panel'+str(panel)]['vertlinewidths'].append(4)
                datafiles['panel'+str(panel)]['vertlinelengths'].append(0.7)
    if datafiles == {}: continue
    for i in range(0,xpanels):
        for j in range(0,ypanels):
            panelnum = i+j*xpanels+1
            ax = fig.add_subplot(ypanels,xpanels,panelnum)
            if 'panel'+str(panelnum) not in datafiles.keys(): continue
            ax.set_xscale(datafiles['panel'+str(panelnum)]['xscale']) 
            ax.set_yscale(datafiles['panel'+str(panelnum)]['yscale']) 
            if datafiles['panel'+str(panelnum)]['ylim'] != []:
                ydown = datafiles['panel'+str(panelnum)]['ylim'][0]
                yup   = datafiles['panel'+str(panelnum)]['ylim'][1]
                ax.set_ylim(ydown,yup)
            else:
                ydown,yup = ax.get_ylim()
            for datafilenum in range(len(datafiles['panel'+str(panelnum)]['lines'].keys())): 
                x = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['x']
                y = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['y']
                linestyle = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['linestyles']
                legend = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['legends']
                drawsty = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['drawsty']
                ax.plot(x,y,linestyle,label=legend,drawstyle=drawsty)
#                if 'peakvelocity' in datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)].keys():
#                    peakvelocity = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['peakvelocity']
#                    peaksnr = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['peaksnr']
#                    ax.vlines(peakvelocity,ydown,yup,linestyle='dotted')
                    #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center',fontsize=20)
            #ax.legend(frameon=False,prop={'size':20},labelspacing=0.1) 
            if j == 0:
                ax.set_title(r'$\rm core'+str(int(corenames[cc]))+r'$')
            ax.text(0.05, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=20)
            ax.text(0.95, 0.9,r'$\rm ('+lletter[j]+')$',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
            if datafiles['panel'+str(panelnum)]['xlim'] != []:
                ax.set_xlim(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1])
                ax.hlines(0,datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],linestyle='dotted')
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            ax.set_xticks(np.arange(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],0.1),minor=True)
            ax.tick_params(which='both',axis='both',direction='in')
            vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
            vertlinexstyles = datafiles['panel'+str(panelnum)]['vertlinestyles']
            vertlinexwidths = datafiles['panel'+str(panelnum)]['vertlinewidths']
            vertlinexlengths = datafiles['panel'+str(panelnum)]['vertlinelengths']
            vertlinexcolors = datafiles['panel'+str(panelnum)]['vertlinecolors']
            for nn,vl in enumerate(vertlinex):
                ax.vlines(vl,ydown,yup*vertlinexlengths[nn],linestyles=vertlinexstyles[nn],colors=vertlinexcolors[nn],linewidths=vertlinexwidths[nn])
            if j != ypanels-1:
                ax.set_yticks(ax.get_yticks()[1:])
                ax.set_xticklabels(ax.get_xlabel(),visible=False)
            else: 
                ax.set_xlabel(xlabel)
            if i != 0:
                ax.set_yticklabels(ax.get_ylabel(),visible=False) 
                ax.set_xticks(ax.get_xticks()[1:]) 
            else: 
                ax.set_ylabel(ylabel)
            
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    plt.close(fig)
    #os.system('open '+pdfname)
    #os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/highres'))
savetxtarr = np.stack((corenames,xw,yw,coremasses,corevelocities[0],coresnr[0],coremom0s[0],coremom1s[0],coremom2s[0],corevelocities[1],coresnr[1],coremom0s[1],coremom1s[1],coremom2s[1],corevelocities[2],coresnr[2],coremom0s[2],coremom1s[2],coremom2s[2],nh3velocities[0],nh3evelocities[0],nh3tkin[0],coregaus_x0[1],coregaus_ex0[1],coregaus_sigma[1],coregaus_esigma[1],coregaus_x0[2],coregaus_ex0[2],coregaus_sigma[2],coregaus_esigma[2]),axis=1)
np.savetxt('convol32_Kirkcores_velocities.txt',savetxtarr,fmt='%3d %10.5f %10.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %5.2f %5.1f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f')
#np.savetxt('single_convol32_Kirkcores_velocities.txt',savetxtarr,fmt='%3d %10.5f %10.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %5.2f %5.1f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f')
#np.savetxt('convol32_Kirkcores_all0p25channel_velocities.txt',savetxtarr,fmt='%3d %10.5f %10.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %5.2f %5.1f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f %7.2f %5.2f')





