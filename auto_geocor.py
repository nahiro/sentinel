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
SCRDIR = os.path.dirname(os.path.abspath(sys.argv[0]))
#REF_DATA_MIN = None
REF_DATA_MIN = 1.0e-5 # for WorldView DN image
RESAMPLING = 'cubic'
RESAMPLING2 = 'near'
MINIMUM_RATIO = 0.9
MINIMUM_NUMBER = 20
REFINE_NUMBER = 10

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog target_georeferenced_image reference_georeferenced_image [options]\n'
'       reference_georeferenced_image is not required if the use_gcps option is given.\n'
'       Both target_georeferenced_image and reference_georeferenced_image are not required if the use_gcps option and the trg_shapefile option are given.\n')
parser.add_option('--scrdir',default=SCRDIR,help='Script directory where find_gcps.py exists (%default)')
parser.add_option('-o','--out_fnam',default=None,help='Output file name (%default)')
parser.add_option('-b','--ref_band',default=None,type='int',help='Reference band# (%default)')
parser.add_option('-B','--trg_band',default=None,type='int',help='Target band# (%default)')
parser.add_option('--ref_multi_band',default=None,type='int',action='append',help='Reference multi-band number (%default)')
parser.add_option('--ref_multi_ratio',default=None,type='float',action='append',help='Reference multi-band ratio (%default)')
parser.add_option('--trg_multi_band',default=None,type='int',action='append',help='Target multi-band number (%default)')
parser.add_option('--trg_multi_ratio',default=None,type='float',action='append',help='Target multi-band ratio (%default)')
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
parser.add_option('--scan_indx_step',default=None,type='int',help='Scan step x index (%default)')
parser.add_option('--scan_indy_step',default=None,type='int',help='Scan step y index (%default)')
parser.add_option('--ref_data_min',default=REF_DATA_MIN,type='float',help='Minimum reference data value (%default)')
parser.add_option('--ref_data_max',default=None,type='float',help='Maximum reference data value (%default)')
parser.add_option('--trg_data_min',default=None,type='float',help='Minimum target data value (%default)')
parser.add_option('--trg_data_max',default=None,type='float',help='Maximum target data value (%default)')
parser.add_option('-r','--rthr',default=None,type='float',help='Threshold of correlation coefficient (%default)')
parser.add_option('-E','--feps',default=None,type='float',help='Step length for curve_fit (%default)')
parser.add_option('--img_fnam',default=None,help='Image file name (%default)')
parser.add_option('-g','--use_gcps',default=None,help='GCP file name to use (%default)')
parser.add_option('-G','--save_gcps',default=None,help='GCP file name to save (%default)')
parser.add_option('--trg_shapefile',default=None,help='Target shapefile (%default)')
parser.add_option('-e','--trg_epsg',default=None,help='Target EPSG (guessed from target data)')
parser.add_option('-n','--npoly',default=None,type='int',help='Order of polynomial used for warping between 1 and 3 (selected based on the number of GCPs)')
parser.add_option('-R','--resampling',default=RESAMPLING,help='Resampling method (%default)')
parser.add_option('--resampling2',default=RESAMPLING2,help='Another resampling method (%default)')
parser.add_option('--resampling2_band',default=None,type='int',action='append',help='Target band# for another resampling method (%default)')
parser.add_option('--minimum_number',default=MINIMUM_NUMBER,type='int',help='Minimum number of GCPs to perform geometric correction (%default)')
parser.add_option('--refine_gcps',default=None,type='float',help='Tolerance to refine GCPs for polynomial interpolation (%default)')
parser.add_option('--minimum_gcps',default=None,type='int',help='Minimum number of GCPs to be left after refine_gcps (available number - discard_number or available number x minimum_ratio)')
parser.add_option('--minimum_ratio',default=MINIMUM_RATIO,type='float',help='Minimum ratio of GCPs to be left after refine_gcps (%default)')
parser.add_option('--discard_number',default=None,type='int',help='Maximum number of GCPs to be discarded by refine_gcps (%default)')
parser.add_option('--refine_number',default=REFINE_NUMBER,type='int',help='Minimum number of GCPs to perform refine_gcps (%default)')
parser.add_option('--tr',default=None,type='float',help='Output resolution in output georeferenced units (%default)')
parser.add_option('--tps',default=False,action='store_true',help='Use thin plate spline transformer (%default)')
parser.add_option('--exp',default=False,action='store_true',help='Output in exp format (%default)')
parser.add_option('--long',default=False,action='store_true',help='Output in long format (%default)')
parser.add_option('-u','--use_edge',default=False,action='store_true',help='Use GCPs near the edge of the correction range (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.use_gcps is None: # Both trg and ref images are required
    if len(args) < 2:
        parser.print_help()
        sys.exit(0)
    trg_fnam = args[0]
    ref_fnam = args[1]
elif opts.trg_shapefile is None: # Only trg image is required
    if len(args) < 1:
        parser.print_help()
        sys.exit(0)
    trg_fnam = args[0]
elif len(args) >= 1: # trg image is optional
    trg_fnam = args[0]
else: # No trg image
    trg_fnam = None

if trg_fnam is not None:
    trg_bnam = os.path.splitext(os.path.basename(trg_fnam))[0]
    tmp_fnam = trg_bnam+'_tmp.tif'
    tmp_xnam = tmp_fnam+'.aux.xml'
    if opts.out_fnam is None:
        out_fnam = trg_bnam+'_geocor.tif'
    else:
        out_fnam = opts.out_fnam
    if opts.trg_epsg is None:
        ds = gdal.Open(trg_fnam)
        prj = ds.GetProjection()
        srs = osr.SpatialReference(wkt=prj)
        opts.trg_epsg = srs.GetAttrValue('AUTHORITY',1)
        ds = None # close dataset

if opts.use_gcps is not None:
    fnam = opts.use_gcps
else:
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'find_gcps.py')
    command += ' '+trg_fnam
    command += ' '+ref_fnam
    command += ' -v'
    if opts.ref_band is not None:
        command += ' --ref_band {}'.format(opts.ref_band)
    if opts.trg_band is not None:
        command += ' --trg_band {}'.format(opts.trg_band)
    if opts.ref_multi_band is not None:
        for band in opts.ref_multi_band:
            command += ' --ref_multi_band {}'.format(band)
    if opts.ref_multi_ratio is not None:
        for ratio in opts.ref_multi_ratio:
            command += ' --ref_multi_ratio {}'.format(ratio)
    if opts.trg_multi_band is not None:
        for band in opts.trg_multi_band:
            command += ' --trg_multi_band {}'.format(band)
    if opts.trg_multi_ratio is not None:
        for ratio in opts.trg_multi_ratio:
            command += ' --trg_multi_ratio {}'.format(ratio)
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
    if opts.scan_indx_step is not None:
        command += ' --scan_indx_step {}'.format(opts.scan_indx_step)
    if opts.scan_indy_step is not None:
        command += ' --scan_indy_step {}'.format(opts.scan_indy_step)
    if opts.ref_data_min is not None:
        command += ' --ref_data_min {}'.format(opts.ref_data_min)
    if opts.ref_data_max is not None:
        command += ' --ref_data_max {}'.format(opts.ref_data_max)
    if opts.trg_data_min is not None:
        command += ' --trg_data_min {}'.format(opts.trg_data_min)
    if opts.trg_data_max is not None:
        command += ' --trg_data_max {}'.format(opts.trg_data_max)
    if opts.rthr is not None:
        command += ' --rthr {}'.format(opts.rthr)
    if opts.feps is not None:
        command += ' --feps {}'.format(opts.feps)
    if opts.img_fnam is not None:
        command += ' --img_fnam {}'.format(opts.img_fnam)
    if opts.exp:
        command += ' --exp'
    if opts.long:
        command += ' --long'
    if opts.use_edge:
        command += ' --use_edge'
    if opts.debug:
        command += ' --debug'
    out = check_output(command,shell=True).decode()
    fnam = StringIO(out)
    if opts.save_gcps is not None:
        with open(opts.save_gcps,'w') as fp:
            fp.write(out)
try:
    xi,yi,xp,yp,dx,dy,r = np.loadtxt(fnam,usecols=(0,1,2,3,4,5,6),unpack=True)
    if xi.size < opts.minimum_number:
        raise ValueError('Error, not enough GCP points.')
except Exception:
    sys.exit()

if trg_fnam is not None:
    command = 'gdal_translate'
    for i,j,x,y in zip(xi,yi,xp,yp):
        command += ' -gcp {} {} {} {}'.format(i,j,x,y)
    command += ' '+trg_fnam
    command += ' '+tmp_fnam
    call(command,shell=True)

    command = 'gdalwarp'
    command += ' -t_srs EPSG:{}'.format(opts.trg_epsg)
    command += ' -overwrite'
    if opts.tr is not None:
        command += ' -tr {} {}'.format(opts.tr,opts.tr)
        command += ' -tap'
    if opts.tps:
        command += ' -tps'
    elif opts.npoly is not None:
        command += ' -order {}'.format(opts.npoly)
    if opts.refine_gcps is not None and xi.size >= opts.refine_number:
        if opts.minimum_gcps is None:
            if opts.discard_number is not None:
                opts.minimum_gcps = xi.size-opts.discard_number
            else:
                opts.minimum_gcps = int(opts.minimum_ratio*xi.size+0.5)
        command += ' -refine_gcps {} {}'.format(opts.refine_gcps,opts.minimum_gcps)
    if opts.resampling2_band is not None:
        tmp2_fnam = trg_bnam+'_tmp2.tif'
        tmp2_xnam = tmp2_fnam+'.aux.xml'
        out1_fnam = trg_bnam+'_geocor1.tif'
        out2_fnam = trg_bnam+'_geocor2.tif'
        command1 = command + ' -r {}'.format(opts.resampling)
        command1 += ' '+tmp_fnam
        command1 += ' '+out1_fnam
        command2 = command + ' -r {}'.format(opts.resampling2)
        command2 += ' '+tmp2_fnam
        command2 += ' '+out2_fnam
        command = 'gdal_translate'
        for band_index in opts.resampling2_band:
            command += ' -b {}'.format(band_index+1)
        command += ' '+tmp_fnam
        command += ' '+tmp2_fnam
        call(command,shell=True)
        call(command1,shell=True)
        call(command2,shell=True)
        ds1 = gdal.Open(out1_fnam)
        ds2 = gdal.Open(out2_fnam)
        drv = gdal.GetDriverByName('GTiff')
        ds = drv.CreateCopy(out_fnam,ds1,strict=0)
        for i,band_index in enumerate(opts.resampling2_band):
            band = ds2.GetRasterBand(i+1)
            ds.GetRasterBand(band_index+1).WriteArray(band.ReadAsArray())
        ds.FlushCache()
        ds = None
        ds1 = None
        ds2 = None
        if os.path.exists(tmp2_fnam):
            os.remove(tmp2_fnam)
        if os.path.exists(tmp2_xnam):
            os.remove(tmp2_xnam)
        if os.path.exists(out1_fnam):
            os.remove(out1_fnam)
        if os.path.exists(out2_fnam):
            os.remove(out2_fnam)
    else:
        command += ' -r {}'.format(opts.resampling)
        command += ' '+tmp_fnam
        command += ' '+out_fnam
        call(command,shell=True)
    if os.path.exists(tmp_fnam):
        os.remove(tmp_fnam)
    if os.path.exists(tmp_xnam):
        os.remove(tmp_xnam)

if opts.trg_shapefile is not None:
    if opts.out_fnam is None:
        out_shapefile = os.path.splitext(os.path.basename(opts.trg_shapefile))[0]+'_geocor.shp'
    else:
        out_shapefile = opts.out_fnam
    command = 'ogr2ogr'
    command += ' -t_srs EPSG:{}'.format(opts.trg_epsg)
    command += ' -overwrite'
    if opts.tps:
        command += ' -tps'
    elif opts.npoly is not None:
        command += ' -order {}'.format(opts.npoly)
    for i,j,x,y in zip(xi,yi,xp,yp):
        command += ' -gcp {} {} {} {}'.format(i,j,x,y)
    command += ' -f "ESRI Shapefile"'
    command += ' '+out_shapefile
    command += ' '+opts.trg_shapefile
    call(command,shell=True)
