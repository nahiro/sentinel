#!/usr/bin/env python
import os
import psutil
mem_size = int(psutil.virtual_memory().available*0.8e-6)
os.environ['_JAVA_OPTIONS'] = '-Xmx{}m'.format(mem_size)     # Save memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx{}m'.format(mem_size))   # Save memory for JAVA
import sys
import re
from snappy import Product,ProductIO,ProductUtils,GPF,HashMap,WKTReader,jpy
from optparse import OptionParser,IndentedHelpFormatter

# Default values

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-b','--band',default=None,type='int',action='append',help='Band# (%default)')
parser.add_option('-G','--geotiff',default=False,action='store_true',help='GeoTiff mode (%default)')
parser.set_usage('Usage: %prog input_fnam [options]')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
input_fnam = args[0]
safe_flag = False
m = re.search('^[^_]+_[^_]+_([^_]+)_.*.zip$',os.path.basename(input_fnam))
if not m:
    m = re.search('^[^_]+_[^_]+_([^_]+)_.*.SAFE$',os.path.basename(input_fnam))
    if not m:
        raise ValueError('Error in file name >>> '+input_fnam)
    safe_flag = True
bnam = os.path.splitext(os.path.basename(input_fnam))[0]
if opts.geotiff:
    output_fnam = '{}.tif'.format(bnam)
else:
    output_fnam = '{}.dim'.format(bnam)
if os.path.exists(output_fnam):
    sys.exit()

# Get snappy Operators
GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
# Read original product
if safe_flag:
    data = ProductIO.readProduct(os.path.join(input_fnam,'MTD_MSIL2A.xml'))
else:
    data = ProductIO.readProduct(input_fnam)
if opts.band is not None:
    band_list = list(data.getBandNames())
    bands = []
    for band in opts.band:
        bands.append(band_list[band])
    params = HashMap()
    params.put('sourceBands',','.join(bands))
    data_tmp = GPF.createProduct('BandSelect',params,data)
    data = data_tmp
if opts.geotiff:
    ProductIO.writeProduct(data,output_fnam,'GeoTiff')
else:
    ProductIO.writeProduct(data,output_fnam,'BEAM-DIMAP')
