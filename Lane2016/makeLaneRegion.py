import numpy as np

corenames, xw, yw, peak850, flux850, cmaj, cmin, cpa = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores.txt',usecols=(0,1,2,3,4,5,6,7),unpack=True,dtype='string')

ff = open('LanecoreRegion.reg','w')
ff.write('# Region file format: DS9 version 4.1\nglobal color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=1 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')

for cc,vv in enumerate(corenames):
    ff.write('ellipse('+xw[cc]+','+yw[cc]+','+str(float(cmaj[cc])/2.)+'",'+str(float(cmin[cc])/2.)+'",'+cpa[cc]+') # text={'+vv+'}\n')

ff.close()

