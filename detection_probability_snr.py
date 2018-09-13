import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import pyfits
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from matplotlib import rc
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

def getbindata(ffalma,ffmirex,mirexfac,contrms,snr,printbins=0):
    hdu_alma = pyfits.open(ffalma)[0] 
    hdu_mirex = pyfits.open(ffmirex)[0]
    alma_data = hdu_alma.data
    mirex_data = hdu_mirex.data*mirexfac # in g cm-2
    print alma_data.shape
    print mirex_data.shape
    
    minSigma = 0
    maxSigma = 0.78
    number_of_bins = 40.
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
        if printbins == 1:
            print 'bb,prob,err',bb,prob,err
        errors.append(err) 
        bbins.append(bb) 
    return bbins,probabilities,errors

p=plt.figure(figsize=(7,6))
plt.subplots_adjust(top=0.94,bottom=0.13,left=0.13,right=0.96)
ax=p.add_subplot(111)
#plt.xlim(0,150)
## linear scale
if scale == 'lin':
    ffalma = 'OrionKLellipse_Lane_on_Stefan_header_CASA.fits'
    ffmirex = 'OrionKLellipse_mask_emap_Orion_A_bw1.0.fits'
    mirexfac = 8.93/214.
    #contrms = 4.69e-4 # from Kirk paper table
    contrms = 10.e-3
    snr = 142.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr)
    ksx1 = np.copy(x)
    ksy1 = np.copy(y)
    plt.errorbar(x,y,yerr=yerror,fmt='b.',ecolor='b',markersize=10,capthick=1.5,zorder=2,label=r'Orion A 15 K')
    # use 20 K for Kirk17 cores
    ffalma = 'OrionKLellipse_Lane_on_Stefan_header_CASA.fits'
    ffmirex = 'OrionKLellipse_mask_emap_Orion_A_bw1.0.fits'
    mirexfac = 8.93/214.
    contrms = 10.e-3
    snr = 142./0.636
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr)
    plt.errorbar(x,y,yerr=yerror,fmt='c.',ecolor='c',capthick=1.5,zorder=2,label=r'Orion A 20 K')
    # alma in K18
    ffalma = 'rebin1p2_pbcor2_uvtaper_briggs_IRDC_C_calibrated_final_cont_2015_image.fits'
    ffmirex = 'wadiao_mirex_on_alma_header_uvtaper.fits'
    mirexfac = 1.
    contrms = 2.e-4
    snr = 3.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr)
    plt.errorbar(x,y,yerr=yerror,fmt='k.',ecolor='k',markersize=10,capthick=1.5,zorder=3,alpha=0.5,label=r'IRDC G28.37')
    # smooth4p8 alma from K18
    ffalma = 'convol4p8_rebin1p2_pbcor2_uvtaper_briggs_IRDC_C_calibrated_final_cont_2015_image.fits'
    ffmirex = 'smooth4p8_wadiao_mirex_on_alma_header_uvtaper.fits'
    mirexfac = 1.
    contrms = 2.e-4*(4.8/2.)**2
    snr = 3.
    x,y,yerror = getbindata(ffalma,ffmirex,mirexfac,contrms,snr)
    ksx2 = np.copy(x)
    ksy2 = np.copy(y)
    plt.errorbar(x,y,yerr=yerror,fmt='r.',ecolor='r',capthick=1.5,zorder=3,alpha=0.5,label=r'IRDC G28.37 smoothed')
    plt.ylim(-0.1,1.1)
    plt.xlim(0,0.8)
    ax.legend(labelspacing=0.1,loc=4,fontsize=12)
    plt.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which='both',axis='both',direction='in')
    ax.text(0.05, 0.95,r'(a)',horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
    pdfname = 'dpf_Orion_G28.pdf'
## logscale
if scale == 'log':
    fitindleft = int(sys.argv[4])
    fitindright = int(sys.argv[5])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.04,1)
    plt.ylim(0.001,2)
    plt.scatter(bbins, probabilities, s = 3, c = 'k', edgecolors = 'face')
    ax.text(0.05, 0.95,panel+r'~~~$\rm SNR\geq'+str(snr)+'$',horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
    xdata = np.array(bbins[fitindleft:fitindright])
    ydata = np.array(probabilities[fitindleft:fitindright])
    yerr = np.array(errors[fitindleft:fitindright])
    pinit = [1.0, 1.0]
    index,indexErr,amp,ampErr = powerlawfit(xdata,ydata,yerr,pinit)
    print 'power',index,'amplitude',amp
    print 'power error',indexErr,'amplitude error',ampErr
    powerlaw = lambda x, amp, index: amp * (x**index) 
    ax.plot(xdata, powerlaw(xdata,amp,index),'b--', label=r'$\rm b='+to_precision(index,4)+r'\pm '+to_precision(indexErr,3)+r'$' '\n' r'$\rm a='+to_precision(amp,4)+r'\pm '+to_precision(ampErr,3)+r'$')
    ax.legend(frameon=False,labelspacing=0.1,loc=3)
    pdfname = 'dpf_snr'+str(snr)+'_log.pdf'
###  
#print len(mirex_data),mirex_data
plt.ylabel(r'$\rm detection~probability$')
plt.xlabel(r'$\rm \Sigma~(g~cm^{-2})$')
os.system('rm '+pdfname)
plt.savefig(pdfname)
os.system('open '+pdfname)
plt.close(p)
os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE/'))

from scipy import stats

ksind1 = (ksx1>=0.08) & (ksx1<=0.28)
ksind2 = (ksx2>=0.08) & (ksx2<=0.28)
newksx1 = ksx1[ksind1]
newksy1 = ksy1[ksind1]
newksx2 = ksx2[ksind2]
newksy2 = ksy2[ksind2]

print newksx1 
print newksy1 
print newksx2 
print newksy2 

ks_d, ks_p = stats.ks_2samp(newksy1,newksy2)
print 'ks_d, ks_p', ks_d, ks_p

