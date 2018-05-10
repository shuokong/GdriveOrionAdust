import numpy as np
import pyfits
import math

ppi = 3.14159265359
cc = 2.99792458e10 # cm/s
hh = 6.62607e-27 # erg*s
kb = 1.38065e-16 # erg/K
ttbg = 2.73 # K

def pdbi(band,tsys,ton,bwidth): # returns rms noise for PdBI.
    #band is the band number, tsys is the system temperature in unit K. ton is the on-source time per configuration in seconds. bwidth is the band width in unit Hz.
    jpk=[22,29,35,45] # conversion factor from K to Jy, in unit Jy/K
    eta=[0.9,0.85,0.8,0.7] # additional efficiency factor due to atmospheric phase noise
    band=band-1 # jpk index is from 0 to 3, so band 2 is jpk[1]
    na=6 # number of antenna
    nc=1 # number of configuration
    npol=2 # number of polarization, use dual mode here #
    return jpk[band]*tsys/eta[band]/math.sqrt(na*(na-1)*nc*ton*bwidth*npol)


def primarybeam(frequency,Diameter):
    """calculate primary beam of antenna in arcsec, given Diameter in m and frequency in Hz"""
    return 1.22*cc/frequency/1.e2/Diameter/ppi*180.*3600.

def bpresp(howmanyFWHM): 
    return math.exp(-0.5*(howmanyFWHM*2.35)**2)

def mbckj(mainbeamFWHM,freq): 
    freq = freq * 1.e9
    return 2.*1.38e7*(1.133*(mainbeamFWHM/3600./180.*ppi)**2.)/(cc/freq)**2.

def kjelli(major,minor,freq): 
    """main beam conversion factor Jy/K for ellipse, freq in GHz, major and minor in arcsec"""
    freq = freq * 1.e9
    return 2.*1.38e7*(1.133*(major/3600./180.*ppi)*(minor/3600./180.*ppi))/(cc/freq)**2.

def mbcjk(mainbeamFWHM,freq):
    """main beam conversion factor K/Jy, freq in GHz, mainbeamFWHM in arcsec"""
    freq = freq * 1.e9
    return 1./(2.*1.38e7*(1.133*(mainbeamFWHM/3600./180.*ppi)**2.)/(cc/freq)**2.)

def jkelli(major,minor,freq):
    """main beam conversion factor K/Jy for ellipse, freq in GHz, major and minor in arcsec"""
    freq = freq * 1.e9
    return 1./(2.*1.38e7*(1.133*(major/3600./180.*ppi)*(minor/3600./180.*ppi))/(cc/freq)**2.)

def almapb(freq): # return ALMA primary beam FWHM (arcsec) at given frequency (Hz)
    D = 12 # m
    llambda = cc / freq
    return 1.17 * llambda / D / 100. / ppi * 180. * 60. * 60.

def pbscale(freq,beam,newfreq):
    """give frequency (GHz) and beamsize (arcsec) and new frequency, calculate new beamsize(arcsec)"""
    return float(freq)/float(newfreq)*float(beam)

def beamsolidangle(mainbeamFWHM):
    """give beam FWHM in arcsec, return Gaussian beam solid angle"""
    return 1.133*(mainbeamFWHM/3600./180.*ppi)**2.

def cubestats(cubefits):
    """get cube rms, max, bmaj, bmin, bpa. rms is the median value in emission free channel (assuming at least 20 emission-free channels). assume cube has stokes, vel, y, x dimensions"""
    hdulist = pyfits.open(cubefits)
    prihdu = hdulist[0]
    data = prihdu.data
    rmscube = prihdu.data[0,:,:,:]
    planermss = np.nanstd(rmscube,axis=(1,2))
    planerms = sorted(planermss)[:20]
    rms = np.nanmedian(planerms)
    cubemax = np.nanmax(data)
    header = prihdu.header
    bmaj = header['BMAJ']*3600.
    bmin = header['BMIN']*3600.
    bpa = header['BPA']
    hdulist.close()
    print 'cube rms =',rms,'cube max =',cubemax,'cube bmaj,bmin,bpa =',bmaj,'arcsec',bmin,'arcsec',bpa,'deg'
    return rms,cubemax,bmaj,bmin,bpa

