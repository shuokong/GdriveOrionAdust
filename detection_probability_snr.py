import sys
import numpy as np
import matplotlib.pyplot as plt
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

print sys.argv
snr = int(sys.argv[1])
scale = sys.argv[2]
panel = sys.argv[3]
print snr, scale


#hdu_alma = pyfits.open('Lane_on_Stefan_header.fits')[0] 
#hdu_mirex = pyfits.open('emap_Orion_A_bw1.0.fits')[1]
hdu_alma = pyfits.open('Lane_on_Stutz_header.fits')[0] 
hdu_mirex = pyfits.open('herschelAmelia/OrionA_all_spire250_nh_mask_corr_apex.fits')[0]

alma_data = hdu_alma.data
alma_header = hdu_alma.header
#mirex_data = hdu_mirex.data*1.67e22/4.27e23 # in g cm-2
mirex_data = hdu_mirex.data/4.27e23 # in g cm-2
mirex_header = hdu_mirex.header

print alma_data.shape
print mirex_data.shape

#masked_mirex = np.copy(mirex_data)
#masked_alma = np.copy(alma_data)
#ymax,xmax = mirex_data.shape
#circlex = 146
#circley = 179
#circler = 25.7
#for y in range(ymax):
#    for x in range(xmax):
#        if ((y-circley)**2 + (x-circlex)**2)**0.5 > circler:
#            masked_mirex[y,x] = -100
#            masked_alma[y,x] = -100
#        else:
#            if masked_mirex[y,x] <= 0.58 or masked_mirex[y,x] >= 0.60:
#                masked_mirex[y,x] = -100
#                masked_alma[y,x] = -100
#            else:
#                if masked_alma[y,x] <= 3.:
#                    masked_alma[y,x] = -100
#
#fits.writeto('masked_alma.fits',masked_alma,header=alma_header,clobber=True)
#fits.writeto('masked_mirex.fits',masked_mirex,header=mirex_header,clobber=True)
#sys.exit()

minSigma = 0
#minSigma = 0.58
maxSigma = 0.78
#maxSigma = 0.64
number_of_bins = 40.
#number_of_bins = 3.
deltaSigma = (maxSigma - minSigma) / number_of_bins
bins = np.arange(minSigma+deltaSigma/2.,maxSigma,deltaSigma)
#print bins
#sys.exit()
contrms = 3.e-4
thres = snr*contrms # in unit of sigma, above which considered detection
probabilities = []
errors = []
bbins = []
for bb in bins:
    mask = (mirex_data<bb+deltaSigma/2.) & (mirex_data>=bb-deltaSigma/2.)
    #print bb,bb-deltaSigma,bb+deltaSigma
    #for y in range(ymax):
    #    for x in range(xmax):
    #        if mask[y,x] == True:
    #            print y,x
    temporary_alma = alma_data[mask] # alma values where mirex in [bb-deltaSigma/2.,bb+deltaSigma/2.)
    temporary1_alma = temporary_alma[~np.isnan(temporary_alma)] # remove nan pixels from alma data
    if len(temporary1_alma) == 0: continue
    temporary2_alma = temporary1_alma[(temporary1_alma>=thres)] # alma snr >= thres is considered detection, this can be changed
    prob = float(len(temporary2_alma))/float(len(temporary1_alma)) # probability = number of alma detection / number of mirex pixels in the bin
    #sys.exit()
    probabilities.append(prob)
    err = (prob*(1.-prob)/float(len(temporary1_alma)))**0.5 
    print len(temporary_alma),len(temporary2_alma),len(temporary1_alma),err
    errors.append(err) 
    bbins.append(bb) 

print probabilities
p=plt.figure(figsize=(7,6))
plt.subplots_adjust(top=0.94,bottom=0.13,left=0.13,right=0.96)
ax=p.add_subplot(111)
#plt.xlim(0,150)
## linear scale
if scale == 'lin':
    x,y = (bbins,probabilities)
    yerror = errors
    plt.errorbar(x,y,yerr=yerror,fmt='k.',ecolor='k',capthick=1.5)
    plt.ylim(-0.1,1.1)
    plt.xlim(0,0.8)
    ax.text(0.5, 0.95,panel+r'~~~$\rm SNR\geq'+str(snr)+'$',horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
    pdfname = 'dpf_snr'+str(snr)+'.pdf'
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
#os.system('cp '+pdfname+os.path.expandvars(' ${DROPATH}/cloudc1'))

