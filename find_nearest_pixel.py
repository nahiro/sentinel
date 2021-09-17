#!/usr/bin/env python
import sys
import gdal
import numpy as np
from optparse import OptionParser,IndentedHelpFormatter

n_nearest = 120

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]

ds = gdal.Open(input_fnam)
data = ds.ReadAsArray()
if ds.RasterCount < 2:
    data_shape = data.shape
else:
    data_shape = data[0].shape
data_trans = ds.GetGeoTransform()
indy,indx = np.indices(data_shape)
xg = data_trans[0]+(indx+0.5)*data_trans[1]+(indy+0.5)*data_trans[2]
yg = data_trans[3]+(indx+0.5)*data_trans[4]+(indy+0.5)*data_trans[5]
ds = None

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
