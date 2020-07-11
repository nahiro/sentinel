#!/usr/bin/env python
import sys
import numpy as np

n_nearest = 120

xstp = 10.0
ystp = -10.0
xmin,xmax,ymin,ymax = (743800.0,756800.0,9236000.0,9251800.0)
xg,yg = np.meshgrid(np.arange(xmin,xmax+0.1*xstp,xstp),np.arange(ymax,ymin-0.1*ystp,ystp))
ngrd = xg.size
nx = xg.shape[1]
ny = xg.shape[0]
sid = np.arange(ngrd).reshape(xg.shape)

sid_0 = []
for nn in range(1,n_nearest+1):
    exec('sid_{} = []'.format(nn))
    exec('leng_{} = []'.format(nn))
for n in range(ngrd):
    i,j = np.unravel_index(n,xg.shape)
    x1 = max(0,j-7)
    x2 = min(nx,x1+15)
    x1 = x2-15
    y1 = max(0,i-7)
    y2 = min(ny,y1+15)
    y1 = y2-15
    l2 = (np.square(xg[y1:y2,x1:x2]-xg[i,j])+np.square(yg[y1:y2,x1:x2]-yg[i,j])).flatten()
    inds = sid[y1:y2,x1:x2].flatten()
    indx = np.argsort(l2)
    if inds[indx[0]] != n:
        raise ValueError('Error, indx[0]={}, n={}'.format(indx[0],n))
    sid_0.append(n)
    for nn in range(1,n_nearest+1):
        exec('sid_{}.append(inds[indx[{}]])'.format(nn,nn))
        exec('leng_{}.append(np.sqrt(l2[indx[{}]]))'.format(nn,nn))
    #break
sid_0 = np.array(sid_0)
for nn in range(1,n_nearest+1):
    exec('sid_{} = np.array(sid_{})'.format(nn,nn))
    exec('leng_{} = np.array(leng_{})'.format(nn,nn))
command = 'np.savez("find_nearest.npz",'
command += 'sid_0=sid_0,'
for nn in range(1,n_nearest+1):
    command += 'sid_{}=sid_{},leng_{}=leng_{},'.format(nn,nn,nn,nn)
command += ')'
exec(command)
