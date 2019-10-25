#!/usr/bin/env python
import os
import sys
import numpy as np
import gdal
import osr
try:
    from io import StringIO
except Exception:
    from StringIO import StringIO
from subprocess import check_output,call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
#REF_DATA_MIN = None
REF_DATA_MIN = 0.1 # for WorldView DN image
RESAMPLING = 'cubic'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog reference_georeferenced_image target_georeferenced_image [options]')
parser.add_option('-b','--ref_band',default=None,type='int',help='Reference band# (%default)')
parser.add_option('-B','--trg_band',default=None,type='int',help='Target band# (%default)')
parser.add_option('-x','--trg_indx_start',default=None,type='int',help='Target start x index (0)')
parser.add_option('-X','--trg_indx_stop',default=None,type='int',help='Target stop x index (target width)')
parser.add_option('-s','--trg_indx_step',default=None,type='int',help='Target step x index (half of subset_width)')
parser.add_option('-y','--trg_indy_start',default=None,type='int',help='Target start y index (0)')
parser.add_option('-Y','--trg_indy_stop',default=None,type='int',help='Target stop y index (target height)')
parser.add_option('-S','--trg_indy_step',default=None,type='int',help='Target step y index (half of subset_height)')
parser.add_option('-W','--subset_width',default=None,type='int',help='Subset width in target pixel (%default)')
parser.add_option('-H','--subset_height',default=None,type='int',help='Subset height in target pixel (%default)')
parser.add_option('--shift_width',default=None,type='int',help='Max shift width in target pixel (%default)')
parser.add_option('--shift_height',default=None,type='int',help='Max shift height in target pixel (%default)')
parser.add_option('--margin_width',default=None,type='int',help='Margin width in target pixel (%default)')
parser.add_option('--margin_height',default=None,type='int',help='Margin height in target pixel (%default)')
parser.add_option('--ref_data_min',default=REF_DATA_MIN,type='float',help='Minimum reference data value (%default)')
parser.add_option('--ref_data_max',default=None,type='float',help='Maximum reference data value (%default)')
parser.add_option('-g','--use_gcps',default=None,help='GCP file name to use (%default)')
parser.add_option('-G','--save_gcps',default=None,help='GCP file name to save (%default)')
parser.add_option('-e','--trg_epsg',default=None,help='Target EPSG (guessed from target data)')
parser.add_option('-r','--rthr',default=None,type='float',help='Threshold of correlation coefficient (%default)')
parser.add_option('-E','--feps',default=None,type='float',help='Step length for curve_fit (%default)')
parser.add_option('-n','--npoly',default=None,type='int',help='Order of polynomial used for warping between 1 and 3 (selected based on the number of GCPs)')
parser.add_option('-R','--resampling',default=RESAMPLING,help='Resampling method (%default)')
parser.add_option('--refine_gcps',default=None,type='float',help='Tolerance to refine GCPs for polynomial interpolation (%default)')
parser.add_option('--tps',default=False,action='store_true',help='Use thin plate spline transformer (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
ref_fnam = args[0]
trg_fnam = args[1]
trg_bnam = os.path.splitext(os.path.basename(trg_fnam))[0]
tmp_fnam = trg_bnam+'_tmp.tif'
out_fnam = trg_bnam+'_geocor.tif'

if opts.trg_epsg is None:
    ds = gdal.Open(trg_fnam)
    prj = ds.GetProjection()
    srs = osr.SpatialReference(wkt=prj)
    opts.trg_epsg = srs.GetAttrValue('AUTHORITY',1)
    ds = None # close dataset

if opts.use_gcps is not None:
    fnam = opts.use_gcps
else:
    command = 'find_gcps.py'
    command += ' '+ref_fnam
    command += ' '+trg_fnam
    command += ' -v'
    if opts.ref_band is not None:
        command += ' --ref_band {}'.format(opts.ref_band)
    if opts.trg_band is not None:
        command += ' --trg_band {}'.format(opts.trg_band)
    if opts.trg_indx_start is not None:
        command += ' --trg_indx_start {}'.format(opts.trg_indx_start)
    if opts.trg_indx_stop is not None:
        command += ' --trg_indx_stop {}'.format(opts.trg_indx_stop)
    if opts.trg_indx_step is not None:
        command += ' --trg_indx_step {}'.format(opts.trg_indx_step)
    if opts.trg_indy_start is not None:
        command += ' --trg_indy_start {}'.format(opts.trg_indy_start)
    if opts.trg_indy_stop is not None:
        command += ' --trg_indy_stop {}'.format(opts.trg_indy_stop)
    if opts.trg_indy_step is not None:
        command += ' --trg_indy_step {}'.format(opts.trg_indy_step)
    if opts.subset_width is not None:
        command += ' --subset_width {}'.format(opts.subset_width)
    if opts.subset_height is not None:
        command += ' --subset_height {}'.format(opts.subset_height)
    if opts.shift_width is not None:
        command += ' --shift_width {}'.format(opts.shift_width)
    if opts.shift_height is not None:
        command += ' --shift_height {}'.format(opts.shift_height)
    if opts.margin_width is not None:
        command += ' --margin_width {}'.format(opts.margin_width)
    if opts.margin_height is not None:
        command += ' --margin_height {}'.format(opts.margin_height)
    if opts.ref_data_min is not None:
        command += ' --ref_data_min {}'.format(opts.ref_data_min)
    if opts.ref_data_max is not None:
        command += ' --ref_data_max {}'.format(opts.ref_data_max)
    if opts.rthr is not None:
        command += ' --rthr {}'.format(opts.rthr)
    if opts.feps is not None:
        command += ' --feps {}'.format(opts.feps)
    if opts.debug:
        command += ' --debug'
    out = check_output(command,shell=True).decode()
    fnam = StringIO(out)
    if opts.save_gcps is not None:
        with open(opts.save_gcps,'w') as fp:
            fp.write(out)
xi,yi,xp,yp,dx,dy,r = np.loadtxt(fnam,unpack=True)

command = 'gdal_translate'
for i,j,x,y in zip(xi,yi,xp,yp):
    command += ' -gcp {} {} {} {}'.format(i,j,x,y)
command += ' '+trg_fnam
command += ' '+tmp_fnam
call(command,shell=True)

command = 'gdalwarp'
command += ' -t_srs EPSG:{}'.format(opts.trg_epsg)
command += ' -overwrite'
if opts.tps:
    command += ' -tps'
elif opts.npoly is not None:
    command += ' -order {}'.format(opts.npoly)
command += ' -r {}'.format(opts.resampling)
if opts.refine_gcps is not None:
    command += ' -refine_gcps {}'.format(opts.refine_gcps)
command += ' '+tmp_fnam
command += ' '+out_fnam
call(command,shell=True)
if os.path.exists(tmp_fnam):
    os.remove(tmp_fnam)
