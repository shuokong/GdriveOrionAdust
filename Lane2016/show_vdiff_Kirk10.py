import math
import sys
import os
from scipy.optimize import curve_fit
from scipy.stats import norm
import numpy as np
import telescope_sensitivity as ts
import pyfits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
from astropy.io import fits
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':20,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

def gaus(x,a,x0,sigma):
    return a/sigma*np.exp(-(x-x0)**2/(2*sigma**2))

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

yso = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/GAScores.txt',usecols=(6),unpack=True,dtype='string')
ysoyes = (yso == 'Y')
ysono = (yso == 'N')
#savetxtarr = np.stack((corenames,xw,yw,coremasses,corevelocities[0],coresnr[0],corevelocities[1],coresnr[1],corevelocities[2],coresnr[2],nh3velocities[0]),axis=1)
corenames,xw,yw,coremasses,corevelocities12CO,coresnr12CO,coremom0_12CO,coremom1_12CO,coremom2_12CO,corevelocities13CO,coresnr13CO,coremom0_13CO,coremom1_13CO,coremom2_13CO,corevelocitiesC18O,coresnrC18O,coremom0_C18O,coremom1_C18O,coremom2_C18O,nh3velocities,nh3evelocities,nh3tkin,coregaus_x0_13CO,coregaus_ex0_13CO,coregaus_sigma_13CO,coregaus_esigma_13CO,coregaus_x0_C18O,coregaus_ex0_C18O,coregaus_sigma_C18O,coregaus_esigma_C18O = np.loadtxt('convol32_Kirkcores_velocities.txt',unpack=True)
#print np.stack((corevelocities13CO.T,corevelocitiesC18O.T),axis=1)
coreysoyes = corenames[ysoyes]
coreysono = corenames[ysono]
print len(coreysoyes),len(coreysono)
#sys.exit()

linenames = [r'$\rm ^{12}CO(1$-$0)$',r'$\rm ^{13}CO(1$-$0)$',r'$\rm C^{18}O(1$-$0)$']
xpanels = 1
ypanels = len(linenames)
xpanelwidth = 10
ypanelwidth = 5
vlow = -3
vhigh = 3
mlow = np.nanmin(coremasses)
mhigh = np.nanmax(coremasses)
cols = len(corenames)

diffvelocities = [nh3velocities-coregaus_x0_C18O,nh3velocities-coregaus_x0_13CO,nh3velocities-coregaus_x0_C18O] # first element dummy
#print 'nh3velocities',nh3velocities
removeind = (nh3velocities<=0.)|(nh3velocities>=16.)|np.isnan(nh3velocities)|(coregaus_x0_13CO<=0.)|(coregaus_x0_13CO>=16.)|np.isnan(coregaus_x0_13CO)|(coregaus_x0_C18O<=0.)|(coregaus_x0_C18O>=16.)|np.isnan(coregaus_x0_C18O)
#print 'removeind',removeind
#sys.exit()
datafiles = {}
j=1
print linenames[j]
panel = 2
print 'panel',panel 
coreveldiff = diffvelocities[panel-1][~removeind]
#print 'len(coreveldiff)',len(coreveldiff)
#print 'np.nanmean(coreveldiff)',np.nanmean(coreveldiff)
#print 'np.nanstd(coreveldiff)',np.nanstd(coreveldiff[abs(coreveldiff)<2])
#print 'scipy.stats.norm.fit(coreveldiff)',norm.fit(coreveldiff[abs(coreveldiff)<2])
#print 'min max coreveldiff',min(coreveldiff),max(coreveldiff)
#ss = raw_input()
hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(vlow,vhigh))
print 'bin size',bin_edges[1]-bin_edges[0]
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
perr = np.sqrt(np.diag(pcov))
print 'popt',popt
print 'perr',perr
datafiles['panel'+str(panel-1)] = {'title':linenames[j],'lines':{'1':{'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'all','linestyles':'k-','drawsty':'steps-mid'},'2':{'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'k--','drawsty':'default'}},'xlim':[vlow,vhigh],'ylim':[0,50],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm NH_3}-v_{\rm gauss}~\rm (km~s^{-1})$','ylabel':r'$\rm number$','text':'','vertlines':[-0.31,0.31]}
coreveldiff = diffvelocities[panel-1][(ysoyes)&(~removeind)]
coreveldiffysoyes = diffvelocities[panel-1][(ysoyes)&(~removeind)]
print 'len(coreveldiff)',len(coreveldiff)
hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(vlow,vhigh))
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
perr = np.sqrt(np.diag(pcov))
print 'popt',popt
print 'perr',perr
datafiles['panel'+str(panel-1)]['lines']['3'] = {'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'YSO','linestyles':'b-','drawsty':'steps-mid'}
datafiles['panel'+str(panel-1)]['lines']['4'] = {'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'b--','drawsty':'default'}
coreveldiff = diffvelocities[panel-1][(ysono)&(~removeind)]
coreveldiffysono = diffvelocities[panel-1][(ysono)&(~removeind)]
print 'len(coreveldiff)',len(coreveldiff)
hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(vlow,vhigh))
bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
perr = np.sqrt(np.diag(pcov))
print 'popt',popt
print 'perr',perr
datafiles['panel'+str(panel-1)]['lines']['5'] = {'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'no YSO','linestyles':'y-','drawsty':'steps-mid'}
datafiles['panel'+str(panel-1)]['lines']['6'] = {'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'y--','drawsty':'default'}

fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels))
gs = gridspec.GridSpec(2,1,height_ratios=[1,2])
plt.subplots_adjust(wspace=0.001,hspace=0.001)
pdfname='corespectra/vdiffgaussKirk10.pdf'
i = 0
panelnum = 1
ax = plt.subplot(gs[panelnum-1])
#ax.set_xscale(datafiles['panel'+str(panelnum)]['xscale']) 
#ax.set_yscale(datafiles['panel'+str(panelnum)]['yscale']) 
if datafiles['panel'+str(panelnum)]['ylim'] != []:
    ydown = datafiles['panel'+str(panelnum)]['ylim'][0]
    yup   = datafiles['panel'+str(panelnum)]['ylim'][1]
    ax.set_ylim(ydown,yup)
if datafiles['panel'+str(panelnum)]['xlim'] != []:
    xmin,xmax = datafiles['panel'+str(panelnum)]['xlim']
    ax.set_xlim(xmin,xmax)
    #ax.hlines(0,xmin,xmax,linestyle='dotted')
for datafilenum in range(len(datafiles['panel'+str(panelnum)]['lines'].keys())): 
    x = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['x']
    y = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['y']
    #xx = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['velocity']
    #print 'weird cores',corenames[(xx<-2)|(xx>2)]
    #ax.hist(x,bins='auto',range=(xmin,xmax))
    linestyle = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['linestyles']
    legend = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['legends']
    drawsty = datafiles['panel'+str(panelnum)]['lines'][str(datafilenum+1)]['drawsty']
    ax.plot(x,y,linestyle,label=legend,drawstyle=drawsty)
    #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center',fontsize=20)
ax.legend(frameon=False,prop={'size':20},labelspacing=0.2) 
#if j == 0:
#    ax.set_title('Gauss')
ax.text(0.05, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
#ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=20)
#ax.text(0.95, 0.9,'('+str(cc+1)+lletter[j]+')',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontsize=20)
xlabel = datafiles['panel'+str(panelnum)]['xlabel']
ylabel = datafiles['panel'+str(panelnum)]['ylabel']
vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
for vl in vertlinex:
    ax.vlines(vl,ydown,yup*0.9,linestyles='dotted',colors='k')
ax.set_yticks(ax.get_yticks()[1:])
ax.set_xticklabels(ax.get_xlabel(),visible=False)
if i != 0:
    ax.set_yticklabels(ax.get_ylabel(),visible=False) 
    ax.set_xticks(ax.get_xticks()[1:]) 
else: 
    ax.set_ylabel(ylabel)
#minor_locator = AutoMinorLocator(5)
#ax.xaxis.set_minor_locator(minor_locator)
ax.tick_params(axis='both',direction='in',length=5,which='major',top=True,right=True)
#ax.tick_params(axis='both',direction='in',length=2,which='minor',top=True,right=True)

############################## for pv diagram

raw_ra_list = np.array([83.52329236,83.57316692,83.62030445,83.6670613,83.71382045,83.75252857,83.79620141,83.82111418,83.84564682,83.86251383,83.86210348,83.85289673,83.84292146,83.82604083,83.80992492,83.80935538,83.80014058,83.79092437,83.77326557,83.7632777,83.74561138,83.73638688,83.73637222,83.74557319,83.75323795,83.77088618,83.79621755,83.82001665,83.84612392,83.8783804,83.91332917,83.95142695,83.98952916,84.0308125,84.07400417,84.11910176,84.1642031,84.21046505]) # same as pvmap_orion.py !
raw_dec_list = np.array([-4.876200218,-4.879461026,-4.897820323,-4.915411506,-4.93299871,-4.965857804,-4.99067018,-5.033833723,-5.076613716,-5.123974925,-5.173573767,-5.222459592,-5.271344718,-5.318700697,-5.36452826,-5.412955261,-5.461988054,-5.510106789,-5.555931851,-5.605577113,-5.6519109,-5.701555547,-5.752730272,-5.801156183,-5.84927762,-5.897401317,-5.940943148,-5.986011267,-6.029551074,-6.068507208,-6.103567778,-6.135131552,-6.166692778,-6.195094167,-6.219704167,-6.240522195,-6.261336664,-6.282571724]) # same as pvmap_orion.py !

xwysoyes = xw[(ysoyes)&(~removeind)]
ywysoyes = yw[(ysoyes)&(~removeind)]
offsetysoyes = []
for nnn,iii in enumerate(xwysoyes):
    xxww = iii
    yyww = ywysoyes[nnn]
    dist = [((xxww-ii)**2+(yyww-raw_dec_list[nn])**2)**0.5 for nn,ii in enumerate(raw_ra_list)]
    offsetysoyes.append(42.-np.argmin(dist)*3.) # raw_ra_list first one at 42 arcmin offset (north), 3 arcmin step pv cut
offsetysoyesarr = np.array(offsetysoyes)

xwysono = xw[(ysono)&(~removeind)]
ywysono = yw[(ysono)&(~removeind)]
offsetysono = []
for nnn,iii in enumerate(xwysono):
    xxww = iii
    yyww = ywysono[nnn]
    dist = [((xxww-ii)**2+(yyww-raw_dec_list[nn])**2)**0.5 for nn,ii in enumerate(raw_ra_list)]
    offsetysono.append(42.-np.argmin(dist)*3.) # raw_ra_list first one at 42 arcmin offset (north), 3 arcmin step pv cut
offsetysonoarr = np.array(offsetysono)

name_cube = 'pv_mask_imfit_13co_pix_2_Tmb_trans_shift.fits'
name_cube_halfmax = 'pv_mask_imfit_13co_pix_2_Tmb_trans_shift_halfmax.fits'

mmap=pyfits.open(name_cube)
cube=mmap[0].data
cdelt1=mmap[0].header['cdelt1'] # velocity
naxis1=mmap[0].header['naxis1']
crval1=mmap[0].header['crval1']
cdelt2=mmap[0].header['cdelt2'] # position
naxis2=mmap[0].header['naxis2']
crval2=mmap[0].header['crval2']
velticks = [-3.,-2.,-1.,0.,1.,2.,3.]
velpix = [(ii-crval1)/cdelt1 for ii in velticks]
veltickslatex = [r"$"+str(ii)+r"$" for ii in velticks]
posticks = [42,21,0,-21,-42,-63]
pospix = [(float(ii)-crval2)/cdelt2 for ii in posticks]
postickslatex = [r"$"+str(ii)+r"$" for ii in posticks]
mmap_halfmax=pyfits.open(name_cube_halfmax)
cube_halfmax=mmap_halfmax[0].data

panelnum = 2
ax2 = plt.subplot(gs[panelnum-1])
ax2.imshow(cube,cmap='gray_r',aspect='auto',extent=[0,naxis1+0.5,naxis2-0.5,-0.5],interpolation='bicubic')
plt.contour(cube_halfmax,levels=[0,1],colors='gray')
#ax2.set_xticklabels(ax2.get_xlabel(),visible=True)
ax2.scatter((coreveldiffysoyes-crval1)/cdelt1,(offsetysoyesarr-crval2)/cdelt2,facecolors='b')
ax2.scatter((coreveldiffysono-crval1)/cdelt1,(offsetysonoarr-crval2)/cdelt2,facecolors='y')
ax2.set_xticks(velpix)
ax2.set_xticklabels(veltickslatex)
ax2.set_yticks(pospix)
ax2.set_yticklabels(postickslatex)
ax2.tick_params(axis='both',direction='in',length=5,which='major',top=True,right=True)
ax2.set_xlabel(xlabel)
ax2.set_ylabel(r'$\rm Offsets~(arcmin)$')
        
os.system('rm '+pdfname)
plt.savefig(pdfname,bbox_inches='tight')
plt.close(fig)
os.system('open '+pdfname)
os.system('cp '+pdfname+os.path.expandvars(' ~/GoogleDrive/imagesSFE/'))


