#!/usr/bin/env python
import os
import sys
import numpy as np
from scipy.interpolate import bisplrep,bisplev
from optparse import OptionParser,IndentedHelpFormatter

# Default values
XTHR = 4.0
YTHR = 4.0
TRG_INDX_STEP = 50
TRG_INDY_STEP = 50
SMOOTH = 1.0e4
DATDIR = os.curdir

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('--xthr',default=XTHR,type='float',help='Max difference in X (%default)')
parser.add_option('--ythr',default=YTHR,type='float',help='Max difference in Y (%default)')
parser.add_option('-s','--trg_indx_step',default=TRG_INDX_STEP,type='int',help='Target step x index (%default)')
parser.add_option('-S','--trg_indy_step',default=TRG_INDY_STEP,type='int',help='Target step y index (%default)')
parser.add_option('--smooth_x',default=SMOOTH,type='float',help='Smoothing factor for X from 0 to 1 (%default)')
parser.add_option('--smooth_y',default=None,type='float',help='Smoothing factor for Y from 0 to 1 (%default)')
parser.add_option('-D','--datdir',default=DATDIR,help='Output data directory (%default)')
parser.add_option('-r','--replace',default=False,action='store_true',help='Replace mode (%default)')
parser.add_option('-e','--exp',default=False,action='store_true',help='Output in exp format (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args
if opts.smooth_y is None:
    opts.smooth_y = opts.smooth_x

for fnam in fnams:
    bnam = os.path.basename(fnam)
    onam = os.path.join(opts.datdir,os.path.splitext(bnam)[0]+'_selected.dat')
    xc,yc,xp,yp,xd,yd,rr = np.loadtxt(fnam,unpack=True)
    xc_uniq = np.unique(xc)
    yc_uniq = np.unique(yc)
    xi = xc.astype(np.int32)
    yi = yc.astype(np.int32)
    xi_uniq = np.unique(xi)
    yi_uniq = np.unique(yi)
    xi_uniq_step = np.diff(xi_uniq).min()
    yi_uniq_step = np.diff(yi_uniq).min()
    xc_offset = xc_uniq[0]-xi_uniq[0]
    yc_offset = yc_uniq[0]-yi_uniq[0]
    if not np.all(xc_uniq[1:]-xi_uniq[1:] == xc_offset):
        raise ValueError('Error, different xc offset >>> {}'.format(fnam))
    if not np.all(yc_uniq[1:]-yi_uniq[1:] == yc_offset):
        raise ValueError('Error, different yc offset >>> {}'.format(fnam))
    xmin = xi.min()
    xmax = xi.max()
    xstp = opts.trg_indx_step
    if xstp != xi_uniq_step:
        sys.stderr.write('Warning, xstp={}, xi_uniq_step={} >>> {}\n'.format(xstp,xi_uniq_step,fnam))
    ymin = yi.min()
    ymax = yi.max()
    ystp = opts.trg_indy_step
    if ystp != yi_uniq_step:
        sys.stderr.write('Warning, ystp={}, yi_uniq_step={} >>> {}\n'.format(ystp,yi_uniq_step,fnam))
    xg,yg = np.meshgrid(np.arange(xmin,xmax+1,xstp),np.arange(ymin,ymax+1,ystp))
    indx = (xi-xmin)//xstp
    indy = (yi-ymin)//ystp
    xc_grid = np.full(xg.shape,np.nan)
    yc_grid = np.full(yg.shape,np.nan)
    xp_grid = np.full(xg.shape,np.nan)
    yp_grid = np.full(yg.shape,np.nan)
    xd_grid = np.full(xg.shape,np.nan)
    yd_grid = np.full(yg.shape,np.nan)
    rr_grid = np.full(xg.shape,np.nan)
    xc_grid[indy,indx] = xc
    yc_grid[indy,indx] = yc
    xp_grid[indy,indx] = xp
    yp_grid[indy,indx] = yp
    xd_grid[indy,indx] = xd
    yd_grid[indy,indx] = yd
    rr_grid[indy,indx] = rr

    xs_grid = np.full(xg.shape,np.nan)
    ys_grid = np.full(yg.shape,np.nan)
    xs_grid[:] = bisplev(np.arange(xg.shape[0]),np.arange(xg.shape[1]),bisplrep(indy,indx,xd,s=opts.smooth_x))
    ys_grid[:] = bisplev(np.arange(yg.shape[0]),np.arange(yg.shape[1]),bisplrep(indy,indx,yd,s=opts.smooth_y))
    cnd = np.isnan(xd_grid)
    xs_grid[cnd] = np.nan
    cnd = np.isnan(yd_grid)
    ys_grid[cnd] = np.nan

    xr_grid = xd_grid-xs_grid
    yr_grid = yd_grid-ys_grid

    with open(onam,'w') as fp:
        for ix,iy in zip(indx,indy):
            if (np.abs(xr_grid[iy,ix]) > opts.xthr) or (np.abs(yr_grid[iy,ix]) > opts.ythr):
                continue
            if opts.replace:
                xp_out = xp_grid[iy,ix]-xd_grid[iy,ix]+xs_grid[iy,ix]
                yp_out = yp_grid[iy,ix]-yd_grid[iy,ix]+ys_grid[iy,ix]
                xd_out = xs_grid[iy,ix]
                yd_out = ys_grid[iy,ix]
            else:
                xp_out = xp_grid[iy,ix]
                yp_out = yp_grid[iy,ix]
                xd_out = xd_grid[iy,ix]
                yd_out = yd_grid[iy,ix]
            if opts.exp:
                line = '{:8.1f} {:8.1f} {:15.8e} {:15.8e} {:15.8e} {:15.8e} {:8.3f}\n'.format(xc_grid[iy,ix],yc_grid[iy,ix],xp_out,yp_out,xd_out,yd_out,rr_grid[iy,ix])
            else:
                line = '{:8.1f} {:8.1f} {:8.2f} {:8.2f} {:6.2f} {:6.2f} {:8.3f}\n'.format(xc_grid[iy,ix],yc_grid[iy,ix],xp_out,yp_out,xd_out,yd_out,rr_grid[iy,ix])
            fp.write(line)
            if opts.verbose:
                sys.stdout.write(line)
