import aplpy
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

mom0 = 1
histogram = 0

hdu1 = fits.open('carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
mask_hdu = fits.open('mask_carmanro_OrionA_all_spire250_nh_mask_corr_apex.fits')[0]
xw, yw = np.loadtxt('../Lane2016/GAScores.txt',usecols=(1,2),unpack=True)

if mom0 == 1:
    xcenter = 83.79424081
    ycenter = -5.578964122
    wid = 0.9290123
    hei = 1.5358419
    xpanels = 1
    ypanels = 1
    fig=plt.figure(figsize=(3*xpanels*1.1*(wid/(wid+hei))*10.,3*ypanels/1.1*(hei/(wid+hei))*10.))
    ff = aplpy.FITSFigure(hdu1, figure=fig)
    ff.recenter(xcenter,ycenter,width=wid,height=hei) 
    ff.set_theme('publication')
    #ff.set_system_latex(True)
    maxcolor = np.nanmax(hdu1.data)
    #maxcolor = 100
    #ff.show_colorscale(cmap='gist_heat', vmin=mincolor, vmax=maxcolor, stretch='log')
    ff.show_colorscale(cmap='gray_r', vmin=1.8e20, vmax=1.4e24, stretch='log')
    ff.show_regions('northfil.reg')
    ff.show_regions('southfil.reg')
    ff.show_regions('olay1.reg')
    ff.show_contour(mask_hdu, levels=1, colors='yellow', linewidths=1, linestyles='dashed', zorder=5)
    ff.show_circles(xw,yw,radius=20./3600.,edgecolor='blue')
    ff.add_colorbar() 
    ff.colorbar.set_font(size=12)
    ff.colorbar.set_pad(0.5)
    ff.colorbar.set_axis_label_text('cm$^{-2}$')
    ff.colorbar.set_font(size=12)
    ff.set_tick_labels_font(size=12)
    ff.set_axis_labels_font(size=12)
    ff.add_scalebar(0.286,corner='bottom right',pad=10) # degree for 2pc at 400 pc
    ff.scalebar.set_label('2 pc')
    ff.scalebar.set_font_size(12)
    #beamx = 83.41442439
    #beamy = -7.022846568
    #bmaj = hdu1.header['BMAJ']
    #bmin = hdu1.header['BMIN']
    #beamangle = hdu1.header['BPA']
    #ff.show_ellipses(beamx,beamy,bmaj,bmin,angle=beamangle-90,facecolor='black',edgecolor='black')
    #ff.add_label(beamx+1.0,beamy+2.0,'0th-moment C$^{18}$O(1-0)',size=12,weight='bold')
    pdfname = 'filregions.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesSFE'))

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
    
    plt.xlabel(r'$I_{\rm ^{12}CO}~\rm K~km~s^{-1}$')
    plt.ylabel('probability density')
    #plt.xlim(0,20)
    ax.xaxis.set_tick_params(top='on',labeltop='on',direction='in')
    ax.yaxis.set_tick_params(direction='in')
    #plt.grid(True)
    pdfname = 'mom0_c18o_hist.pdf'
    os.system('rm '+pdfname)
    plt.savefig(pdfname,bbox_inches='tight')
    os.system('open '+pdfname)
    os.system('cp '+pdfname+os.path.expandvars(' /Users/shuokong/GoogleDrive/imagesCARMAOrion/'))

