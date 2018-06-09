import math
import sys
import os
import numpy as np
import telescope_sensitivity as ts
import pyfits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':12,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

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
corenames, xw, yw, coremasses = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/GAScores.txt',usecols=(0,1,2,3),unpack=True)
worldcoord = np.stack((xw,yw,np.zeros_like(xw),np.zeros_like(xw)),axis=1)
#worldcoord = np.stack((xw,yw,np.zeros_like(xw)),axis=1)

cellsize = 2. # voxel size in arcsec
#lines = ['/Users/shuokong/GoogleDrive/12co/products/mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/13co/products/mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/c18o/products/mask_imfit_c18o_pix_2_Tmb.fits']
lines = ['/Users/shuokong/GoogleDrive/12co/products/convol32_mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/13co/products/convol32_mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/c18o/products/convol32_mask_imfit_c18o_pix_2_Tmb.fits'] # convol32 to match GAS survey
#lines = ['/Users/shuokong/GoogleDrive/Alyssa/nostokes_mask_imfit_12co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/Alyssa/regrid_12co_specsmooth_0p25_mask_imfit_13co_pix_2_Tmb.fits','/Users/shuokong/GoogleDrive/Alyssa/regrid_12co_specsmooth_0p25_mask_imfit_c18o_pix_2_Tmb.fits']
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
    #scidata = hdulist[0].data[:,:,:] # usually index order: v, dec, ra
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
    pixcoord = w.all_world2pix(worldcoord,1) # FITS standard uses 1
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
print 'finish getting nh3 metadata'

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
coresnr = [[],[],[]]
for cc in range(cols):
    fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels/1.1))
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    #pdfname='corespectra/averspec_Kirk_core'+str(cc+1)+'.pdf'
    pdfname='corespectra/convol32_averspec_Kirk_core'+str(cc+1)+'.pdf'
    #pdfname='corespectra/averspec_all0p25channel_core'+str(cc+1)+'.pdf'
    datafiles = {}
    #print 'x,y',x,y
    ccnh3xx = int(nh3xx[cc])
    ccnh3yy = int(nh3yy[cc])
    try:
        ccvlsr = [nh3data[ccnh3yy,ccnh3xx]]
    except:
        ccvlsr = [-100]
    print 'ccnh3xx,ccnh3yy',ccnh3xx,ccnh3yy,'ccvlsr',ccvlsr
    nh3velocities[0].append(ccvlsr[0])
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
                continue
            rawspectrum = np.nanmean(temp,axis=(1))
            crval3, cdelt3, crpix3, n1 = line3rdaxis[panel-1]
            rawvelocity = np.array([(crval3+cdelt3*((ii+1)-crpix3))/1.e3 for ii in range(n1)]) # should not use i, it will confuse with the for i above
            subvel = (rawvelocity > vlow) & (rawvelocity < vhigh) 
            rawintens = rawspectrum[subvel]
            velocity = rawvelocity[subvel]
            coresavetxtarr = np.stack((velocity,rawintens),axis=1)
            np.savetxt('corespectra/datapoints/'+molnames[j]+'/convol32_Kirkcore'+str(int(corenames[cc]))+'_'+molnames[j]+'.txt',coresavetxtarr,fmt='%7.2f %10.2e')
            coremom0 = mom0(cdelt3/1.e3,rawintens)
            coremom0s[panel-1].append(coremom0)
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
            datafiles['panel'+str(panel)] = {'title':linenames[j],'lines':{'1':{'x':velocity,'y':rawintens,'peakvelocity':velocity[peakind],'peaksnr':snr,'legends':'data','linestyles':'k-','drawsty':'steps-mid'}},'xlim':[vlow,vhigh],'ylim':[ymin-(ymax-ymin)/10.,ymax+(ymax-ymin)/10.],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm LSR}~\rm (km~s^{-1})$','ylabel':r'$T_{\rm mb}~\rm (K)$','text':'','vertlines':[ccvlsr[0],coremom1],'vertlinestyles':[ccvlsrstyle[0],'dotted']}
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
                peakvelocity = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['peakvelocity']
                peaksnr = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['peaksnr']
                linestyle = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['linestyles']
                legend = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['legends']
                drawsty = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['drawsty']
                ax.plot(x,y,linestyle,label=legend,drawstyle=drawsty)
                ax.vlines(peakvelocity,ydown,yup,linestyle='dashed')
                ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center',fontsize=12)
            #ax.legend(frameon=False,prop={'size':14},labelspacing=0.1) 
            if j == 0:
                ax.set_title('core'+str(int(corenames[cc])))
            ax.text(0.05, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            ax.text(0.95, 0.9,'('+str(cc+1)+lletter[j]+')',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            if datafiles['panel'+str(panelnum)]['xlim'] != []:
                ax.set_xlim(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1])
                ax.hlines(0,datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],linestyle='dotted')
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            ax.set_xticks(np.arange(datafiles['panel'+str(panelnum)]['xlim'][0],datafiles['panel'+str(panelnum)]['xlim'][1],0.1),minor=True)
            vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
            vertlinexstyles = datafiles['panel'+str(panelnum)]['vertlinestyles']
            for nn,vl in enumerate(vertlinex):
                ax.vlines(vl,ydown,yup*0.6,linestyles=vertlinexstyles[nn],colors='b')
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
savetxtarr = np.stack((corenames,xw,yw,coremasses,corevelocities[0],coresnr[0],coremom0s[0],coremom1s[0],coremom2s[0],corevelocities[1],coresnr[1],coremom0s[1],coremom1s[1],coremom2s[1],corevelocities[2],coresnr[2],coremom0s[2],coremom1s[2],coremom2s[2],nh3velocities[0]),axis=1)
#np.savetxt('Kirkcores_velocities.txt',savetxtarr,fmt='%3d %10.5f %10.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f')
np.savetxt('convol32_Kirkcores_velocities.txt',savetxtarr,fmt='%3d %10.5f %10.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f')
#np.savetxt('Lanecores_all0p25channel_peak_velocities.txt',savetxtarr,fmt='%3d %15.5f %15.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f')




