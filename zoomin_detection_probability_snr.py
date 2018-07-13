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

def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

def powerlawfit(xdata,ydata,yerr,pinit): # xdata,ydata,yerr n-element arrays, pinit two-element list

    ##########
    # Fitting the data -- Least Squares Method
    ##########

    # xdata = np.linspace(1.1, 10.1, num_points)
    # ydata = powerlaw(xdata, 10.0, -2.0)     # simulated perfect data
    # yerr = 0.2 * ydata
    
    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    # Note that the `logyerr` term here is ignoring a constant prefactor.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)
    #
    
    logx = np.log10(xdata)
    logy = np.log10(ydata)
    logyerr = yerr / ydata
    
    # define our (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    
    # pinit = [1.0, -1.0]
    out = optimize.leastsq(errfunc, pinit, args=(logx, logy, logyerr), full_output=1)
    
    pfinal = out[0]
    covar = out[1]
    print pfinal
    print covar
    
    index = pfinal[1]
    amp = 10.0**pfinal[0]
    
    indexErr = np.sqrt( covar[1][1] )
    ampErr = np.sqrt( covar[0][0] ) * amp 
     
    return index,indexErr,amp,ampErr 

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
    return bbins,probabilities,errors

def coremass(S850,kappa=0.012,Td=15.,D=400.):
    """S850 in Jy, kappa in cm2 g-1, Td in K, D in pc"""
    return 1.3*S850*(kappa/0.012)*(np.exp(17./Td)-1.)*(D/450.)**2

corenames, xw, yw, peak850, flux850, cmaj, cmin, cpa = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores_degree.txt',usecols=(0,1,2,3,4,5,6,7),unpack=True)
print 'max(cmaj)',max(cmaj)
print 'max(cmin)',max(cmin)
print 'sum(cmaj>30)',sum(cmaj>30)
print 'sum(cmin>30)',sum(cmin>30)
print 'sum(cmaj>60)',sum(cmaj>60)
print 'sum(cmin>60)',sum(cmin>60)
c = SkyCoord(xw, yw, "icrs", unit="deg")

mass850 = coremass(flux850)
print 'type(mass850)',type(mass850)
print 'len(mass850) in Msun',len(mass850)
print 'np.nansum(mass850) in Msun',np.nansum(mass850) # need to scale down by a factor of 0.8 because Lane uses 450 pc

#ffalma = 'Lane_on_Stefan_header.fits'
ffalma = 'Lane_on_Stefan_header_CASA.fits'
ffmirex = 'mask_emap_Orion_A_bw1.0.fits'
hdu_jcmt = pyfits.open('Lane2016/OrionA_850_auto_mos_clip.fits')[0]
print 'np.nansum(hdu_jcmt.data[hdu_jcmt.data>4.69e-4*5.])',np.nansum(hdu_jcmt.data[hdu_jcmt.data>4.69e-4*5.])
print 'jcmt SNR>5 flux/np.nansum(flux850)',np.nansum(hdu_jcmt.data[hdu_jcmt.data>4.69e-4*5.])/np.nansum(flux850)
hdu_nicest = pyfits.open(ffmirex)[0]
print 'np.nansum(hdu_nicest.data)',np.nansum(hdu_nicest.data)
print 'np.nansum(hdu_nicest.data) to mass in Msun',np.nansum(hdu_nicest.data)*1.67e22/4.27e23*(30.*400.*1.5e13)**2/2.e33
print 'continuum mass fraction',np.nansum(mass850)/(np.nansum(hdu_nicest.data)*1.67e22/4.27e23*(30.*400.*1.5e13)**2/2.e33)
print 'jcmt SNR>5 mass fraction',np.nansum(hdu_jcmt.data[hdu_jcmt.data>4.69e-4*5.])/np.nansum(flux850)*np.nansum(mass850)/(np.nansum(hdu_nicest.data)*1.67e22/4.27e23*(30.*400.*1.5e13)**2/2.e33)
#sys.exit()
header = hdu_nicest.header
w = wcs.WCS(header)
scidata = hdu_nicest.data
pixcoordx,pixcoordy = c.galactic.to_pixel(w)
xx = pixcoordx.astype(int)
yy = pixcoordy.astype(int)
coreak = scidata[yy,xx]
mirexfac = 9.
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
rawhist, bin_edges = np.histogram(coreav,bins=range(0,61,2),range=(0,60))
hist = rawhist.astype(float)
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
print 'bincenter',bincenter
rawav_hist, av_bin_edges = np.histogram(mapav,bins=range(0,61,2),range=(0,60))
av_hist = rawav_hist.astype(float)
dpf_hist = hist/av_hist
print 'dpf_hist',dpf_hist
print 'type(hist)',type(hist)
print 'type(av_hist)',type(av_hist)
#sys.exit()
#contrms = 4.69e-4 # from Kirk paper table
contrms = 10.e-3

p=plt.figure(figsize=(7,12))
plt.subplots_adjust(top=0.98,bottom=0.1,left=0.12,right=0.96,hspace=0.1)
pdfname = 'dpf_OrionA.pdf'

ax=p.add_subplot(211)
avmin,avmax = (0,60)
nnbins = (avmax-avmin)/2
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
    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(avmin,avmax)
    ax.legend(frameon=False,labelspacing=0.5,loc=2,handletextpad=0.5,fontsize=12)
    ax.text(0.98, 0.95,'(a)',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
    ax.set_ylabel(r'$\rm detection~probability$')
###  
ax=p.add_subplot(212)
avmin,avmax = (0,20)
nnbins = (avmax-avmin)/2
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
    ax.set_ylim(-0.1,0.6)
    ax.set_xlim(avmin,avmax)
    ax.legend(frameon=False,labelspacing=0.5,loc=2,handletextpad=0.5,fontsize=12)
    ax.text(0.98, 0.95,'(b)',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
    ax.set_ylabel(r'$\rm detection~probability$')
    ax.set_xlabel(r'$\rm A_V~(mag)$')

os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
os.system('open '+pdfname)
plt.close(p)
#os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/cloudc1'))

