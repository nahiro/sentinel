#!/usr/bin/env python
import sys
import numpy as np

n_nearest = 120

sid,xc,yc,ndat,leng,area = np.loadtxt('get_center.dat',unpack=True)
sid = (sid+0.1).astype(np.int64)
ndat = (ndat+0.1).astype(np.int64)
sid_indx = np.arange(sid.size)

sid_0 = []
for nn in range(1,n_nearest+1):
    exec('sid_{} = []'.format(nn))
    exec('leng_{} = []'.format(nn))
for n in range(sid.size):
    l2 = np.square(xc-xc[n])+np.square(yc-yc[n])
    indx = np.where((l2 < 1.0e-4) & (sid_indx != n))[0]
    if len(indx) != 0:
        sys.stderr.write('Warning, n={:6d}'.format(n))
        for i in indx:
            sys.stderr.write(', l2[{:6d}]={}'.format(i,l2[i]))
        sys.stderr.write('\n')
    indx = np.argsort(l2)
    #if sid[indx[0]] != n:
    #    #raise ValueError('Error, indx[0]={}, n={}'.format(indx[0],n))
    #    sys.stderr.write('Warning, indx[0]={}, n={}\n'.format(indx[0],n))
    sid_0.append(n)
    for nn in range(1,n_nearest+1):
        exec('sid_{}.append(sid[indx[{}]])'.format(nn,nn))
        exec('leng_{}.append(np.sqrt(l2[indx[{}]]))'.format(nn,nn))
    #break
sid_0 = np.array(sid_0)
for nn in range(1,n_nearest+1):
    exec('sid_{} = np.array(sid_{})'.format(nn,nn))
    exec('leng_{} = np.array(leng_{})'.format(nn,nn))
command = 'np.savez("find_nearest_field.npz",'
command += 'sid_0=sid_0,'
for nn in range(1,n_nearest+1):
    command += 'sid_{}=sid_{},leng_{}=leng_{},'.format(nn,nn,nn,nn)
command += ')'
exec(command)
