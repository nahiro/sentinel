#!/usr/bin/env python
import os
import psutil
mem_size = int(psutil.virtual_memory().available*0.8e-6)
os.environ['_JAVA_OPTIONS'] = '-Xmx{}m'.format(mem_size)     # Set memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx{}m'.format(mem_size))   # Set memory for JAVA
import sys
import re
from snappy import Product,ProductIO,ProductUtils,GPF,HashMap,WKTReader,jpy
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
ORIGIN_X = 743805.0 # pixel center
ORIGIN_Y = 9236005.0 # pixel center
XY_STEP = 10.0
EPSG = 32748 # UTM zone 48S

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog input_fnam [options]')
parser.add_option('-g','--gamma0',default=False,action='store_true',help='Output gamma0 instead of sigma0 (%default)')
parser.add_option('--skip_orbit',default=False,action='store_true',help='Do not apply orbit file (%default)')
parser.add_option('--speckle',default=False,action='store_true',help='Apply speckle filter (%default)')
parser.add_option('--iangle',default=False,action='store_true',help='Output incidence angle (%default)')
parser.add_option('-x','--origin_x',default=ORIGIN_X,type='float',help='X origin of standard grid (%default)')
parser.add_option('-y','--origin_y',default=ORIGIN_Y,type='float',help='Y origin of standard grid (%default)')
parser.add_option('-s','--xy_step',default=XY_STEP,type='float',help='XY step of standard grid (%default)')
parser.add_option('-E','--epsg',default=EPSG,help='Output EPSG (%default)')
parser.add_option('-S','--std_grid',default=False,action='store_true',help='Use standard grid (%default)')
parser.add_option('-T','--tiff',default=False,action='store_true',help='GeoTiff mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]
m = re.search('^[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*.zip$',os.path.basename(input_fnam))
if not m:
    m = re.search('^[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*.SAFE$',os.path.basename(input_fnam))
    if not m:
        raise ValueError('Error in file name >>> '+input_fnam)
dstr = m.group(1)[:8]
if opts.tiff:
    output_fnam = '{}.tif'.format(dstr)
else:
    output_fnam = '{}.dim'.format(dstr)
if os.path.exists(output_fnam):
    sys.exit()

# Get snappy Operators
GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
# Read original product
data = ProductIO.readProduct(input_fnam)
# Apply orbit file (ApplyOrbitFileOp.java)
if not opts.skip_orbit:
    params = HashMap()
    try:
        data_tmp = GPF.createProduct('Apply-Orbit-File',params,data)
    except Exception:
        sys.stderr.write('Warning, error in applying orbit file >>> '+dstr+'\n')
        data_tmp = data
    data = data_tmp
# Subset (SubsetOp.java)
WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')
wkt = "POLYGON((107.201 -6.910,107.367 -6.910,107.367 -6.760,107.201 -6.760,107.201 -6.910))"
geom = WKTReader().read(wkt)
params = HashMap()
params.put('copyMetadata',True)
params.put('geoRegion',geom)
data_tmp = GPF.createProduct('Subset',params,data)
data = data_tmp
# Calibration (CalibrationOp.java)
params = HashMap()
if opts.gamma0:
    params.put('outputSigmaBand',False)
    params.put('outputBetaBand',True)
else:
    params.put('outputSigmaBand',True)
data_tmp = GPF.createProduct('Calibration',params,data)
data = data_tmp
# Terrain flattening (TerrainFlatteningOp.java)
if opts.gamma0:
    params = HashMap()
    data_tmp = GPF.createProduct('Terrain-Flattening',params,data)
    data = data_tmp
# Speckle reduction (SpeckleFilterOp.java)
if opts.speckle:
    params = HashMap()
    data_tmp = GPF.createProduct('Speckle-Filter',params,data)
    data = data_tmp
# Terrain correction (RangeDopplerGeocodingOp.java)
params = HashMap()
if opts.iangle:
    params.put('saveIncidenceAngleFromEllipsoid',True)
params.put('demName','SRTM 3Sec')
params.put('demResamplingMethod','BILINEAR_INTERPOLATION')
params.put('imgResamplingMethod','BILINEAR_INTERPOLATION')
params.put('pixelSpacingInMeter',opts.xy_step)
params.put('mapProjection','EPSG:{}'.format(opts.epsg))
#params.put('mapProjection','AUTO:42001') # WGS84/AutoUTM
if opts.std_grid:
    params.put('alignToStandardGrid',True)
    params.put('standardGridOriginX',opts.origin_x)
    params.put('standardGridOriginY',opts.origin_y)
data_tmp = GPF.createProduct("Terrain-Correction",params,data)
data = data_tmp
# Convert to dB (LinearTodBOp.java)
params = HashMap()
if opts.iangle:
    band_list = list(data.getBandNames())
    bands = []
    for band in band_list:
        if not re.search('angle',band.lower()):
            bands.append(band)
    params.put('sourceBands',','.join(bands))
data_tmp = GPF.createProduct('linearToFromdB',params,data)
if opts.iangle:
    for band in band_list:
        if re.search('angle',band.lower()):
            ProductUtils.copyBand(band,data,band,data_tmp,True)
data = data_tmp
# Attach bandname (BandSelectOp.java)
band_list = list(data.getBandNames())
bands = []
for band in band_list:
    band_new = band+'_'+dstr
    ProductUtils.copyBand(band,data,band_new,data,True)
    bands.append(band_new)
params = HashMap()
params.put('sourceBands',','.join(bands))
data_tmp = GPF.createProduct('BandSelect',params,data)
data = data_tmp
if opts.tiff:
    ProductIO.writeProduct(data,output_fnam,'GeoTiff')
else:
    ProductIO.writeProduct(data,output_fnam,'BEAM-DIMAP')
