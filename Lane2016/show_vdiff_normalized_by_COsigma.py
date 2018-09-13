import math
import sys
import os
from scipy.optimize import curve_fit
import numpy as np
import telescope_sensitivity as ts
import pyfits
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
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
ratiolow = -3
ratiohigh = 3
mlow = np.nanmin(coremasses)
mhigh = np.nanmax(coremasses)
cols = len(corenames)

usegausvel = 1

if usegausvel == 1:
    diffvelocities = [(nh3velocities-coregaus_x0_C18O)/coregaus_sigma_C18O,(nh3velocities-coregaus_x0_13CO)/coregaus_sigma_13CO,(nh3velocities-coregaus_x0_C18O)/coregaus_sigma_C18O] # first element dummy
    #print 'nh3velocities',nh3velocities
    removeind = (nh3velocities<=0.)|(nh3velocities>=16.)|np.isnan(nh3velocities)|(coregaus_x0_13CO<=0.)|(coregaus_x0_13CO>=16.)|np.isnan(coregaus_x0_13CO)|(coregaus_sigma_13CO<=0.)|(coregaus_sigma_13CO>=16.)|np.isnan(coregaus_sigma_13CO)|(coregaus_x0_C18O<=0.)|(coregaus_x0_C18O>=16.)|np.isnan(coregaus_x0_C18O)|(coregaus_sigma_C18O<=0.)|(coregaus_sigma_C18O>=16.)|np.isnan(coregaus_sigma_C18O)
    #print 'removeind',removeind
    #sys.exit()
    datafiles = {}
    for i in range(0,xpanels):
        for j in range(1,ypanels): # start from 13CO
            print linenames[j]
            panel = i+j*xpanels+1
            print 'panel',panel 
            coreveldiff = diffvelocities[panel-1][~removeind]
            print 'len(coreveldiff)',len(coreveldiff)
            print 'min max coreveldiff',min(coreveldiff),max(coreveldiff)
            hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(ratiolow,ratiohigh))
            bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
            #popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
            #perr = np.sqrt(np.diag(pcov))
            #print 'popt',popt
            #print 'perr',perr
            datafiles['panel'+str(panel-1)] = {'title':linenames[j],
                'lines':{
                         '1':{'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'all','linestyles':'k-','drawsty':'steps-mid'},
                         #'2':{'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'k--','drawsty':'default'},
                         },
                      'xlim':[ratiolow,ratiohigh],'ylim':[0,np.nanmax(hist)*1.1],'xscale':'linear','yscale':'linear','xlabel':r'$(v_{\rm NH_3}-v_{\rm gauss})/\sigma_{\rm gauss}$','ylabel':r'$\rm number$','text':'','vertlines':[-1,1]}
            coreveldiff = diffvelocities[panel-1][(ysoyes)&(~removeind)]
            print 'len(coreveldiff)',len(coreveldiff)
            hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(ratiolow,ratiohigh))
            bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
            #popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
            #perr = np.sqrt(np.diag(pcov))
            #print 'popt',popt
            #print 'perr',perr
            datafiles['panel'+str(panel-1)]['lines']['3'] = {'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'YSO','linestyles':'b-','drawsty':'steps-mid'}
            #datafiles['panel'+str(panel-1)]['lines']['4'] = {'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'b--','drawsty':'default'}
            coreveldiff = diffvelocities[panel-1][(ysono)&(~removeind)]
            print 'len(coreveldiff)',len(coreveldiff)
            hist, bin_edges = np.histogram(coreveldiff,bins='auto',range=(ratiolow,ratiohigh))
            bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
            #popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
            #perr = np.sqrt(np.diag(pcov))
            #print 'popt',popt
            #print 'perr',perr
            datafiles['panel'+str(panel-1)]['lines']['5'] = {'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'no YSO','linestyles':'y-','drawsty':'steps-mid'}
            #datafiles['panel'+str(panel-1)]['lines']['6'] = {'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':'fit '+r'$\sigma$='+'{0:.3f}'.format(popt[2])+r'$\pm$'+'{0:.3f}'.format(perr[2]),'linestyles':'y--','drawsty':'default'}
    
    fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels/1.1))
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    pdfname='corespectra/vdiff_normalized_by_COsigma.pdf'
    ypanels = 2
    for i in range(0,xpanels):
        for j in range(0,ypanels):
            panelnum = i+j*xpanels+1
            ax = fig.add_subplot(ypanels,xpanels,panelnum)
            if 'panel'+str(panelnum) not in datafiles.keys(): continue
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
            for lll in datafiles['panel'+str(panelnum)]['lines'].keys(): 
                x = datafiles['panel'+str(panelnum)]['lines'][lll]['x']
                y = datafiles['panel'+str(panelnum)]['lines'][lll]['y']
                #xx = datafiles['panel'+str(panelnum)]['lines'][lll]['velocity']
                #print 'weird cores',corenames[(xx<-2)|(xx>2)]
                #ax.hist(x,bins='auto',range=(xmin,xmax))
                linestyle = datafiles['panel'+str(panelnum)]['lines'][lll]['linestyles']
                legend = datafiles['panel'+str(panelnum)]['lines'][lll]['legends']
                drawsty = datafiles['panel'+str(panelnum)]['lines'][lll]['drawsty']
                ax.plot(x,y,linestyle,label=legend,drawstyle=drawsty)
                #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center',fontsize=12)
            ax.legend(frameon=False,labelspacing=0.2) 
            if j == 0:
                ax.set_title('Gauss')
            ax.text(0.05, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            #ax.text(0.95, 0.9,'('+str(cc+1)+lletter[j]+')',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
            for vl in vertlinex:
                ax.vlines(vl,ydown,yup*0.9,linestyles='dashed',colors='k')
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
            minor_locator = AutoMinorLocator(5)
            ax.xaxis.set_minor_locator(minor_locator)
            ax.tick_params(axis='both',direction='in',length=5,which='major',top=True,right=True)
            ax.tick_params(axis='both',direction='in',length=2,which='minor',top=True,right=True)
            
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    plt.close(fig)
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' ~/GoogleDrive/imagesSFE/'))


