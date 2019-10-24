#!/usr/bin/env python
import os
import sys
import numpy as np
try:
    from io import StringIO
except Exception:
    from StringIO import StringIO
from subprocess import check_output,call
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog reference_geotiff_file target_geotiff_file [options]')
parser.add_option('-b','--ref_band',default=None,type='int',help='Reference band# (%default)')
parser.add_option('-B','--trg_band',default=None,type='int',help='Target band# (%default)')
parser.add_option('-W','--subset_width',default=None,type='int',help='Subset width in target pixel (%default)')
parser.add_option('-H','--subset_height',default=None,type='int',help='Subset height in target pixel (%default)')
parser.add_option('--shift_width',default=None,type='int',help='Max shift width in target pixel (%default)')
parser.add_option('--shift_height',default=None,type='int',help='Max shift height in target pixel (%default)')
parser.add_option('--margin_width',default=None,type='int',help='Margin width in target pixel (%default)')
parser.add_option('--margin_height',default=None,type='int',help='Margin height in target pixel (%default)')
parser.add_option('-r','--rthr',default=None,type='float',help='Threshold of correlation coefficient (%default)')
parser.add_option('-E','--feps',default=None,type='float',help='Step length for curve_fit (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
ref_fnam = args[0]
trg_fnam = args[1]
tmp_fnam = 'tmp.tif'
out_fnam = os.path.splitext(os.path.basename(trg_fnam))[0]+'_geocor.tif'

command = 'sen2_geocor.py'
command += ' '+ref_fnam
command += ' '+trg_fnam
command += ' -v'
if opts.ref_band is not None:
    command += ' --ref_band {}'.format(opts.ref_band)
if opts.trg_band is not None:
    command += ' --trg_band {}'.format(opts.trg_band)
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
if opts.rthr is not None:
    command += ' --rthr {}'.format(opts.rthr)
if opts.feps is not None:
    command += ' --feps {}'.format(opts.feps)
if opts.debug:
    command += ' --debug'
out = check_output(command,shell=True)

xi,yi,xp,yp,dx,dy,r = np.loadtxt(StringIO(out.decode()),unpack=True)
command = 'gdal_translate'
for i,j,x,y in zip(xi,yi,xp,yp):
    command += ' -gcp {} {} {} {}'.format(i,j,x,y)
command += ' '+trg_fnam
command += ' '+tmp_fnam
call(command,shell=True)

command = 'gdalwarp'
command += ' -t_srs EPSG:32748'
command += ' -r cubic'
command += ' -order 2'
#command += ' -refine_gcps 0.5'
command += ' -overwrite'
#command += ' -tps'
command += ' '+tmp_fnam
command += ' '+out_fnam
call(command,shell=True)
os.remove(tmp_fnam)
