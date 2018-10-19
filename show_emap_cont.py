import aplpy
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

peak = 0
lane = 1
temp = 0
histogram = 0

hdu1 = fits.open('mask_emap_Orion_A_bw1.0.fits')[0]
hdu2 = fits.open('Lane_on_Stefan_header_CASA.fits')[0]
hdu3 = fits.open('../AncillaryData/GBT/OrionA_Tkin_DR1_rebase3_flag.fits')[0]

if peak == 1:
    xcenter=210.8647791
    ycenter=-19.47488042
    wid = 5.4750000
    hei = 1.8083333
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*(wid/(wid+hei))*10.,3*ypanels*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(hdu1, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    #ff.set_system_latex(True)
    maxcolor = np.nanmax(hdu1.data)
    ff.show_colorscale(cmap='gray_r', vmin=0, vmax=maxcolor, stretch='sqrt')
    ff.show_regions('OrionKLellipse.reg')
    #ff.show_contour(hdu2, levels=0.01*np.concatenate((np.array([5,100,200]),np.arange(300,30000,600))), colors='yellow', linewidths=0.5)
    #ff.show_contour(hdu2, levels=[0.01*10.], colors='yellow', linewidths=0.5)
    ff.show_contour(hdu2, levels=[0.01*142.], colors='yellow', linewidths=1.0)
    ff.show_contour(hdu1, levels=[20./9.], colors='blue', linewidths=0.8)
    ff.add_colorbar() 
    ff.colorbar.set_font(size=20)
    ff.colorbar.set_pad(0.2)
    #ff.colorbar.set_axis_label_text(r'$A_K$')
    ff.set_tick_labels_font(size=20)
    ff.set_axis_labels_font(size=20)
    ff.add_scalebar(0.286,corner='bottom',pad=1) # degree for 2pc at 400 pc
    ff.scalebar.set_label('2 pc') 
    ff.scalebar.set_font_size(20) 
    beamx = 213.59927
    beamy = -20.18993
    bmaj = 60./3600.
    bmin = 60./3600.
    beamangle = 0
    ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    #ff.add_label(beamx+1.0,beamy+2.0,'Peak Intensity $^{12}$CO(1-0)',size=20,weight='bold')
    ff.show_lines([np.array([[210,210],[-18.5,-20.5]])],linestyles='dashed',color='black',zorder=5)
    ff.tick_labels.set_xformat('dd')
    ff.tick_labels.set_yformat('dd.d')
    pdfname = 'emap_cont.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE/'))

if lane == 1:
    xcenter=210.8647791
    ycenter=-19.47488042
    wid = 5.4750000
    hei = 1.8083333
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*(wid/(wid+hei))*10.,3*ypanels*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(hdu2, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    #ff.set_system_latex(True)
    maxcolor = np.nanmax(hdu2.data)
    ff.show_colorscale(cmap='gray_r', vmin=0.01*5, vmax=maxcolor, stretch='log')
    ff.set_tick_labels_font(size=20)
    ff.set_axis_labels_font(size=20)
    ff.add_scalebar(0.286,corner='bottom',pad=1) # degree for 2pc at 400 pc
    ff.scalebar.set_label('2 pc') 
    ff.scalebar.set_font_size(20) 
    beamx = 213.59927
    beamy = -20.18993
    bmaj = 60./3600.
    bmin = 60./3600.
    beamangle = 0
    ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    ff.tick_labels.set_xformat('dd')
    ff.tick_labels.set_yformat('dd.d')
    pdfname = 'lanecont.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE/'))

if temp == 1:
    xcenter=83.89552443
    ycenter=-5.590266212
    wid = 0.5778173
    hei = 1.5565690
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*(wid/(wid+hei))*10.,3*ypanels*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(hdu3, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    #ff.set_system_latex(True)
    maxcolor = np.nanmax(hdu1.data)
    ff.show_colorscale(cmap='gray', vmin=5, vmax=30, stretch='linear')
    ff.show_regions('OrionKLellipse.reg')
    #ff.show_contour(hdu2, levels=[0.01*10.], colors='yellow', linewidths=0.5)
    ff.show_contour(hdu2, levels=[0.01*142.], colors='yellow', linewidths=1.0)
    ff.add_colorbar() 
    ff.colorbar.set_font(size=20)
    ff.colorbar.set_pad(0.2)
    #ff.colorbar.set_axis_label_text(r'$A_K$')
    ff.set_tick_labels_font(size=20)
    ff.set_axis_labels_font(size=20)
    ff.add_scalebar(0.286,corner='bottom right',pad=1) # degree for 2pc at 400 pc
    ff.scalebar.set_label('2 pc') 
    ff.scalebar.set_font_size(20) 
    #beamx = 213.59927
    #beamy = -20.18993
    #bmaj = 60./3600.
    #bmin = 60./3600.
    #beamangle = 0
    #ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    #ff.add_label(beamx+1.0,beamy+2.0,'Peak Intensity $^{12}$CO(1-0)',size=20,weight='bold')
    #ff.show_lines([np.array([[210,210],[-18.5,-20.5]])],linestyles='dashed',color='black',zorder=5)
    ff.tick_labels.set_xformat('dd.d')
    ff.tick_labels.set_yformat('dd.d')
    pdfname = 'tkinmap.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE/'))

if histogram == 1:
    from matplotlib import rc
    rc('text', usetex=True)
    font = {'weight':'normal','size':12,'family':'sans-serif','sans-serif':['Helvetica']}
    rc('font', **font)
    
    print hdu1.data.shape
    #sys.exit()
    x = hdu1.data[~np.isnan(hdu1.data)]
    p=plt.figure(figsize=(7,6))
    fig, ax = plt.subplots(1,1)
    # the histogram of the data
    n, bins, patches = plt.hist(x, 100, normed=True, histtype='step', color='k')
    #print n,bins
    
    plt.xlabel(r'$T_{\rm peak,12}~\rm K$')
    plt.ylabel('probability density')
    #plt.xlim(0,20)
    ax.xaxis.set_tick_params(top='on',labeltop='on',direction='in')
    ax.yaxis.set_tick_params(direction='in')
    #plt.grid(True)
    pdfname = 'peak_12co_hist.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMAOrion/'))
