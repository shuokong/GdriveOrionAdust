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
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rc
from astropy.io import fits
import astropy.wcs as wcs
rc('text', usetex=True)
font = {'weight' : 'normal','size':30,'family':'sans-serif','sans-serif':['Helvetica']}
rc('font', **font)

def gaus(x,a,x0,sigma):
    return a/sigma*np.exp(-(x-x0)**2/(2*sigma**2))

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

#savetxtarr = np.stack((corenames,xw,yw,coremasses,corevelocities[0],coresnr[0],corevelocities[1],coresnr[1],corevelocities[2],coresnr[2],nh3velocities[0]),axis=1)
corenames,xw,yw,coremasses,corevelocities12CO,coresnr12CO,coremom0_12CO,coremom1_12CO,coremom2_12CO,corevelocities13CO,coresnr13CO,coremom0_13CO,coremom1_13CO,coremom2_13CO,corevelocitiesC18O,coresnrC18O,coremom0_C18O,coremom1_C18O,coremom2_C18O,nh3velocities,nh3evelocities,nh3tkin,coregaus_x0_13CO,coregaus_ex0_13CO,coregaus_sigma_13CO,coregaus_esigma_13CO,coregaus_x0_C18O,coregaus_ex0_C18O,coregaus_sigma_C18O,coregaus_esigma_C18O = np.loadtxt('convol32_Kirkcores_velocities.txt',unpack=True)
#print np.stack((corevelocities13CO.T,corevelocitiesC18O.T),axis=1)

removeind = (nh3velocities<=0.)|(nh3velocities>=16.)|np.isnan(nh3velocities)|(coregaus_x0_13CO<=0.)|(coregaus_x0_13CO>=16.)|np.isnan(coregaus_x0_13CO)|(coregaus_x0_C18O<=0.)|(coregaus_x0_C18O>=16.)|np.isnan(coregaus_x0_C18O)
#print 'removeind',removeind
#sys.exit()

nhhdu = fits.open('../herschelAmelia/carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
north_maskhdu = fits.open('../herschelAmelia/northmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
south_maskhdu = fits.open('../herschelAmelia/southmask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]

worldcoord2 = np.stack((xw,yw),axis=1) # 
w = wcs.WCS(nhhdu.header)
pixcoord = w.all_world2pix(worldcoord2,1) # FITS standard uses 1
xx = pixcoord[:,0]
yy = pixcoord[:,1]

topcluster_13CO_vdiffgauss = [nh3velocities[nn]-coregaus_x0_13CO[nn] for nn,ii in enumerate(corenames) if north_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1 and nn not in removeind and abs(nh3velocities[nn]-coregaus_x0_13CO[nn]) < 3]
botcluster_13CO_vdiffgauss = [nh3velocities[nn]-coregaus_x0_13CO[nn] for nn,ii in enumerate(corenames) if south_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1 and nn not in removeind and abs(nh3velocities[nn]-coregaus_x0_13CO[nn]) < 3]

topcluster_C18O_vdiffgauss = [nh3velocities[nn]-coregaus_x0_C18O[nn] for nn,ii in enumerate(corenames) if north_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1 and nn not in removeind and abs(nh3velocities[nn]-coregaus_x0_C18O[nn]) < 3]
botcluster_C18O_vdiffgauss = [nh3velocities[nn]-coregaus_x0_C18O[nn] for nn,ii in enumerate(corenames) if south_maskhdu.data[int(yy[nn]),int(xx[nn])] == 1 and nn not in removeind and abs(nh3velocities[nn]-coregaus_x0_C18O[nn]) < 3]

print 'min max topcluster_13CO_vdiffgauss',min(topcluster_13CO_vdiffgauss),max(topcluster_13CO_vdiffgauss)
print 'min max topcluster_C18O_vdiffgauss',min(topcluster_C18O_vdiffgauss),max(topcluster_C18O_vdiffgauss)
print 'min max botcluster_13CO_vdiffgauss',min(botcluster_13CO_vdiffgauss),max(botcluster_13CO_vdiffgauss)
print 'min max botcluster_C18O_vdiffgauss',min(botcluster_C18O_vdiffgauss),max(botcluster_C18O_vdiffgauss)
print 'scipy.stats.norm.fit(topcluster_13CO_vdiffgauss)',norm.fit(topcluster_13CO_vdiffgauss)
print 'scipy.stats.norm.fit(topcluster_C18O_vdiffgauss)',norm.fit(topcluster_C18O_vdiffgauss)
print 'scipy.stats.norm.fit(botcluster_13CO_vdiffgauss)',norm.fit(botcluster_13CO_vdiffgauss)
print 'scipy.stats.norm.fit(botcluster_C18O_vdiffgauss)',norm.fit(botcluster_C18O_vdiffgauss)
#sys.exit()

linenames = [r'$\rm ^{12}CO(1$-$0)$',r'$\rm ^{13}CO(1$-$0)$',r'$\rm C^{18}O(1$-$0)$']
xpanels = 1
ypanels = len(linenames)
xpanelwidth = 10
ypanelwidth = 5
vlow = -3
vhigh = 3

north = 0
south = 1

### north
if north == 1:
    diffvelocities = [topcluster_13CO_vdiffgauss,topcluster_13CO_vdiffgauss,topcluster_C18O_vdiffgauss] # first element dummy
    #print 'nh3velocities',nh3velocities
    datafiles = {}
    for i in range(0,xpanels):
        for j in range(1,ypanels): # start from 13CO
            print linenames[j]
            panel = i+j*xpanels+1
            print 'panel',panel 
            coreveldiff = diffvelocities[panel-1]
            print 'len(coreveldiff)',len(coreveldiff)
            print 'min max coreveldiff',min(coreveldiff),max(coreveldiff)
            hist, bin_edges = np.histogram(coreveldiff,bins=60,range=(vlow,vhigh))
            bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
            popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
            perr = np.sqrt(np.diag(pcov))
            print 'popt',popt
            print 'perr',perr
            datafiles['panel'+str(panel-1)] = {'title':linenames[j],'lines':{'1':{'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'','linestyles':'k-','drawsty':'steps-mid'},'2':{'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':r'$\rm x_0='+'{0:.2f}'.format(popt[1])+r'\pm'+'{0:.2f}'.format(perr[1])+'$\n'+r'$\sigma='+'{0:.2f}'.format(popt[2])+r'\pm'+'{0:.2f}'.format(perr[2])+r'$','linestyles':'b--','drawsty':'default'}},'xlim':[vlow,vhigh],'ylim':[0,30],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm ce}~\rm (km~s^{-1})$','ylabel':r'$\rm number$','text':'','vertlines':[popt[1]]}
    
    fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels/1.1))
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    pdfname='corespectra/vdiffgauss_northcluster.pdf'
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
                #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center')
            ax.legend(frameon=False,prop={'size':25},loc=2,labelspacing=0.2,handletextpad=0.1)
            if j == 0:
                ax.set_title(r'$\rm North$')
            ax.text(0.95, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            #ax.text(0.95, 0.9,'('+str(cc+1)+lletter[j]+')',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
            for vl in vertlinex:
                ax.vlines(vl,ydown,yup*0.9,linestyles='dashed',colors='b')
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

### south
if south == 1:
    diffvelocities = [botcluster_13CO_vdiffgauss,botcluster_13CO_vdiffgauss,botcluster_C18O_vdiffgauss] # first element dummy
    #print 'nh3velocities',nh3velocities
    datafiles = {}
    for i in range(0,xpanels):
        for j in range(1,ypanels): # start from 13CO
            print linenames[j]
            panel = i+j*xpanels+1
            print 'panel',panel 
            coreveldiff = diffvelocities[panel-1]
            print 'len(coreveldiff)',len(coreveldiff)
            print 'min max coreveldiff',min(coreveldiff),max(coreveldiff)
            hist, bin_edges = np.histogram(coreveldiff,bins=60,range=(vlow,vhigh))
            bincenter = (bin_edges[:-1] + bin_edges[1:]) / 2.
            popt,pcov = curve_fit(gaus,bincenter,hist,p0=[1,0,0.5])
            perr = np.sqrt(np.diag(pcov))
            print 'popt',popt
            print 'perr',perr
            datafiles['panel'+str(panel-1)] = {'title':linenames[j],'lines':{'1':{'x':bincenter,'y':hist,'velocity':coreveldiff,'peaksnr':[],'legends':'','linestyles':'k-','drawsty':'steps-mid'},'2':{'x':bincenter,'y':gaus(bincenter,*popt),'velocity':coreveldiff,'peaksnr':[],'legends':r'$\rm x_0='+'{0:.2f}'.format(popt[1])+r'\pm'+'{0:.2f}'.format(perr[1])+'$\n'+r'$\sigma='+'{0:.2f}'.format(popt[2])+r'\pm'+'{0:.2f}'.format(perr[2])+r'$','linestyles':'b--','drawsty':'default'}},'xlim':[vlow,vhigh],'ylim':[0,30],'xscale':'linear','yscale':'linear','xlabel':r'$v_{\rm ce}~\rm (km~s^{-1})$','ylabel':r'$\rm number$','text':'','vertlines':[popt[1]]}
    
    fig=plt.figure(figsize=(xpanelwidth*xpanels*1.1,ypanelwidth*ypanels/1.1))
    plt.subplots_adjust(wspace=0.001,hspace=0.001)
    pdfname='corespectra/vdiffgauss_southcluster.pdf'
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
                #ax.text(peakvelocity+0.8, yup*0.9, '%.1f' % peakvelocity + ',' + '%.1f' % peaksnr,horizontalalignment='left',verticalalignment='center')
            ax.legend(frameon=False,prop={'size':25},loc=2,labelspacing=0.2,handletextpad=0.1)
            if j == 0:
                ax.set_title(r'$\rm South$')
            ax.text(0.95, 0.9,datafiles['panel'+str(panelnum)]['title'],horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
            #ax.text(0.1, 0.9,datafiles['panel'+str(panelnum)]['title']+' '+datafiles['panel'+str(panelnum)]['text'],horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            #ax.text(0.95, 0.9,'('+str(cc+1)+lletter[j]+')',horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontsize=12)
            xlabel = datafiles['panel'+str(panelnum)]['xlabel']
            ylabel = datafiles['panel'+str(panelnum)]['ylabel']
            vertlinex = datafiles['panel'+str(panelnum)]['vertlines']
            for vl in vertlinex:
                ax.vlines(vl,ydown,yup*0.9,linestyles='dashed',colors='b')
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

