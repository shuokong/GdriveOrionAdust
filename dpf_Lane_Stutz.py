import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import pyfits
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from matplotlib import rc
from astropy.coordinates import SkyCoord
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':20,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)
from scipy import optimize 

scale = 'lin'

def getbindata(ffalma,ffmirex,mirexfac,contrms,snr,minAv,maxAv,nbins):
    hdu_alma = pyfits.open(ffalma)[0] 
#    if 'emap' in ffmirex:
#        hdu_mirex = pyfits.open(ffmirex)[1]
#    else:
    hdu_mirex = pyfits.open(ffmirex)[0]
    alma_data = hdu_alma.data
    mirex_data = hdu_mirex.data*mirexfac # in Av
    print alma_data.shape
    print mirex_data.shape
    
    minSigma = minAv
    maxSigma = maxAv
    number_of_bins = nbins
    deltaSigma = (maxSigma - minSigma) / number_of_bins
    bins = np.arange(minSigma+deltaSigma/2.,maxSigma,deltaSigma)
    contrms = contrms 
    thres = snr*contrms # in unit of sigma, above which considered detection
    probabilities = []
    errors = []
    bbins = []
    for bb in bins:
        mask = (mirex_data<bb+deltaSigma/2.) & (mirex_data>=bb-deltaSigma/2.)
        temporary_alma = alma_data[mask] # alma values where mirex in [bb-deltaSigma/2.,bb+deltaSigma/2.)
        temporary1_alma = temporary_alma[~np.isnan(temporary_alma)] # remove nan pixels from alma data
        if len(temporary1_alma) == 0: continue
        temporary2_alma = temporary1_alma[(temporary1_alma>=thres)] # alma snr >= thres is considered detection, this can be changed
        prob = float(len(temporary2_alma))/float(len(temporary1_alma)) # probability = number of alma detection / number of mirex pixels in the bin
        probabilities.append(prob)
        err = (prob*(1.-prob)/float(len(temporary1_alma)))**0.5 
        print len(temporary_alma),len(temporary2_alma),len(temporary1_alma),err
        errors.append(err) 
        bbins.append(bb) 
    return np.array(bbins),np.array(probabilities),np.array(errors)

def coremass(S850,kappa=0.012,Td=15.,D=400.):
    """S850 in Jy, kappa in cm2 g-1, Td in K, D in pc"""
    return 1.3*S850*(kappa/0.012)*(np.exp(17./Td)-1.)*(D/450.)**2

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
                pixellist.append((x,y))
    return pixellist

##
make_Lanecores_on_Stutz_header = 0 #
if make_Lanecores_on_Stutz_header == 1:
    corenames, xw, yw, peak850, flux850, cmaj, cmin, cpa = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores_degree.txt',usecols=(0,1,2,3,4,5,6,7),unpack=True)
    print 'max(cmaj)',max(cmaj)
    print 'max(cmin)',max(cmin)
    print 'sum(cmaj>30)',sum(cmaj>30)
    print 'sum(cmin>30)',sum(cmin>30)
    print 'sum(cmaj>60)',sum(cmaj>60)
    print 'sum(cmin>60)',sum(cmin>60)
    cc = SkyCoord(xw, yw, "icrs", unit="deg")
    
    ffjcmt = 'Lane2016/nofreq_OrionA_850_auto_mos_clip.fits' # only for core pixels
    hdu_jcmt = pyfits.open(ffjcmt)
    header = hdu_jcmt[0].header
    cellsize = abs(header['CDELT1']*header['CDELT2'])**0.5*3600.
    ww = wcs.WCS(header)
    #scidata = hdu_jcmt[0].data
    #scidata *= 0
    pixcoordx,pixcoordy = cc.to_pixel(ww)
    xx = pixcoordx.astype(int)
    yy = pixcoordy.astype(int)
    #scidata[yy,xx] = 1
    ffemap = 'Lane_on_Stutz_header.fits'
    hdu_emap = pyfits.open(ffemap)
    header_emap = hdu_emap[0].header
    ww_emap = wcs.WCS(header_emap)
    emapdata = hdu_emap[0].data
    emapdata *= 0
    for ccc,vvv in enumerate(corenames):
        pixlist = beampixel(cmaj[ccc],cmin[ccc],cpa[ccc],[xx[ccc],yy[ccc]],cellsize)
        pixlistarr = np.array(pixlist)
        sc_jcmt = SkyCoord.from_pixel(pixlistarr[:,0], pixlistarr[:,1],ww)
        pixcoordx,pixcoordy = sc_jcmt.galactic.to_pixel(ww_emap)
        #scidata[pixlistarr[:,1],pixlistarr[:,0]] = 1
        xxx = np.around(pixcoordx).astype(int)
        yyy = np.around(pixcoordy).astype(int)
        emapdata[yyy,xxx] = 1
    #hdu_jcmt.writeto('Lanecores_nofreq_OrionA_850_auto_mos_clip.fits',clobber=True)
    hdu_emap.writeto('Lanecores_on_Stutz_header.fits',clobber=True)
    hdu_jcmt.close()
    hdu_emap.close()
    sys.exit()
###

ffalma = 'Lane_on_Stutz_header.fits'
ffmirex = 'herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits' 

hdu_nicest = pyfits.open(ffmirex)[0]
scidata = hdu_nicest.data

ff_lanecores = 'Lanecores_on_Stutz_header.fits' # from above

hdu_lanecores = pyfits.open(ff_lanecores)[0]
scidata_lanecores = hdu_lanecores.data
lanecores_detection = (scidata_lanecores>0)
coreak = scidata[lanecores_detection]
mirexfac = 1./2.e21
coreav = coreak * mirexfac
mapav = scidata.flatten() * mirexfac
#print 'coreak[:3]',coreak[:3] # checked in ds9, ok
median_coreav = np.nanmedian(coreav)
coreav_low = np.nanmin(coreav)
coreav_high = np.nanmax(coreav)
print 'coreav_low,coreav_high',coreav_low,coreav_high
print 'median_coreav',median_coreav
print 'sum(coreav<0)',sum(coreav<0)
print 'sum(coreav>0)',sum(coreav>0)
bininterval = 2 # mag
rawhist, bin_edges = np.histogram(coreav,bins=range(0,61,bininterval),range=(0,60))
hist = rawhist.astype(float)
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
print 'bincenter',bincenter
rawav_hist, av_bin_edges = np.histogram(mapav,bins=range(0,61,bininterval),range=(0,60))
av_hist = rawav_hist.astype(float)
dpf_hist = hist/av_hist
print 'dpf_hist',dpf_hist
print 'type(hist)',type(hist)
print 'type(av_hist)',type(av_hist)
#sys.exit()
#contrms = 4.69e-4 # from Kirk paper table
contrms = 2.5e-4

p=plt.figure(figsize=(7,12))
plt.subplots_adjust(top=0.98,bottom=0.1,left=0.12,right=0.96,hspace=0.1)
pdfname = 'dpf_OrionA_Stutz.pdf'

ax=p.add_subplot(211)
avmin,avmax = (0,60)
nnbins = (avmax-avmin)/bininterval
## linear scale
if scale == 'lin':
    snr = 3.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    ax.errorbar(x,y,yerr=yerror,fmt='r.',ecolor='r',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    snr = 5.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    ax.errorbar(x,y,yerr=yerror,fmt='g.',ecolor='g',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    snr = 10.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    fr_10_20 = (x>10)&(x<20)
    ppp = np.polyfit(x[fr_10_20],y[fr_10_20],deg=1,w=yerror[fr_10_20])
    print 'fr_10_20 ppp',ppp
    fr_20_30 = (x>20)&(x<30)
    ppp = np.polyfit(x[fr_20_30],y[fr_20_30],deg=1,w=yerror[fr_20_30])
    print 'fr_20_30 ppp',ppp
    ax.errorbar(x,y,yerr=yerror,fmt='b.',ecolor='b',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    ax.errorbar(bincenter,dpf_hist,yerr=(dpf_hist*(1-dpf_hist)/av_hist)**0.5,drawstyle='steps-mid',color='k',capthick=1.5,zorder=2,label=r'$\rm core$')
    ax.vlines(20,-0.1,0.7,linestyles='dashed')
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(avmin,avmax)
    ax.legend(frameon=False,labelspacing=0.5,loc=2,handletextpad=0.5,fontsize=20)
    ax.text(0.98, 0.95,'(a)',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
    ax.set_ylabel(r'$\rm sub$-$\rm mm~detection~probability$')
###  
bininterval = 1 # mag
rawhist, bin_edges = np.histogram(coreav,bins=range(0,61,bininterval),range=(0,60))
hist = rawhist.astype(float)
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
print 'bincenter',bincenter
rawav_hist, av_bin_edges = np.histogram(mapav,bins=range(0,61,bininterval),range=(0,60))
av_hist = rawav_hist.astype(float)
dpf_hist = hist/av_hist
ax=p.add_subplot(212)
avmin,avmax = (0,20)
nnbins = (avmax-avmin)/bininterval
## linear scale
if scale == 'lin':
    snr = 3.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    ax.errorbar(x,y,yerr=yerror,fmt='r.',ecolor='r',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    snr = 5.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    ax.errorbar(x,y,yerr=yerror,fmt='g.',ecolor='g',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    snr = 10.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr,avmin,avmax,nnbins)
    ax.errorbar(x,y,yerr=yerror,fmt='b.',ecolor='b',capthick=1.5,zorder=2,label=r'$\rm pixel~SNR\geq~$'+str(int(snr)))
    ax.errorbar(bincenter,dpf_hist,yerr=(dpf_hist*(1-dpf_hist)/av_hist)**0.5,drawstyle='steps-mid',color='k',capthick=1.5,zorder=2,label=r'$\rm core$')
    ax.set_ylim(-0.05,0.3)
    ax.set_xlim(avmin+1,avmax-5)
    ax.legend(frameon=False,labelspacing=0.5,loc=2,handletextpad=0.5,fontsize=20)
    ax.text(0.98, 0.95,'(b)',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
    ax.set_ylabel(r'$\rm sub$-$\rm mm~detection~probability$')
    ax.set_xlabel(r'$\rm cloud~A_V~(mag)$')

os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
os.system('open '+pdfname)
plt.close(p)
os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE/'))

