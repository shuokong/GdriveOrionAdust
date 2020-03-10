import sys
import os
import numpy as np

lletter = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

########################
corenames, xw, yw, peak850, flux850, cmaj, cmin, cpa, ssize, proto, ind = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores_degree.txt',usecols=(0,1,2,3,4,5,6,7,8,9,10),unpack=True,dtype='string')

regff = open('Lanecores_starless.reg','w')
regff.write('# Region file format: DS9 version 4.1\n')
regff.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
regff.write('fk5\n')
for nn,ii in enumerate(corenames):
    if proto[nn] == 'N':
        regff.write('circle('+xw[nn]+','+yw[nn]+',15.000") # text={'+ii+'}\n') 
regff.close()

########################
corenames, xw, yw, peak850, flux850, cmaj, cmin, cpa, ssize, proto, ind = np.loadtxt('/Users/shuokong/GoogleDrive/OrionAdust/Lane2016/Getsources_cores.txt',usecols=(0,1,2,3,4,5,6,7,8,9,10),unpack=True,dtype='string')
already = ['99','655','282','357','294','610','255','309']

catff = open('Lanecores_starless.cat','w')
catff.write('#  catalog file for Kong Orion A\n')
catff.write('#  Orion A cloud:  (check LSR velocity)\n')
catff.write('\n')
catff.write('# Lane cores:\n')
for nn,ii in enumerate(corenames):
    if proto[nn] == 'N':
        if ii in already: 
            catff.write('#'+xw[nn]+'   '+yw[nn]+'    J2000.0  Lane'+str(ii).zfill(3)+'  8. LSR RAD\n') 
        else:
            catff.write(xw[nn]+'   '+yw[nn]+'    J2000.0  Lane'+str(ii).zfill(3)+'  8. LSR RAD\n') 
catff.write('\n')
catff.write('#  Ref. positions:\n')
catff.write('5:12:41   -7:25:44    J2000.0  Stick_REF_W   8. LSR RAD\n')
catff.write('#5:29:00   -5:25:30    J2000.0  Stick_REF_N   8. LSR RAD\n')
catff.write('#5:29:16   -6:55:21    J2000.0  Stick_REF_S   8. LSR RAD\n')

catff.close()

