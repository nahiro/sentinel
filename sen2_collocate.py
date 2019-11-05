#!/usr/bin/env python
import os
import psutil
mem_size = int(psutil.virtual_memory().available*0.8e-6)
os.environ['_JAVA_OPTIONS'] = '-Xmx{}m'.format(mem_size)     # Save memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx{}m'.format(mem_size))   # Save memory for JAVA
import sys
import re
from datetime import datetime
from snappy import Product,ProductIO,ProductUtils,GPF,HashMap,jpy
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
FNAM = 'collocate_all.tif'
BAND = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-F','--output_fnam',default=FNAM,help='Output tiff file name (%default)')
parser.add_option('-b','--band',default=None,type='int',action='append',help='Band# ({})'.format(BAND))
parser.add_option('-s','--skip_rename_master',default=False,action='store_true',help='Do not rename master bands (%default)')
parser.add_option('-S','--skip_rename_slave',default=False,action='store_true',help='Do not rename slave bands (%default)')
(opts,args) = parser.parse_args()
if len(args) < 2:
    parser.print_help()
    sys.exit(0)
fnams = args[:]
if opts.band is None:
    opts.band = BAND

# Get snappy Operators
GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
# Read original product
data_1 = ProductIO.readProduct(fnams[0])
# Attach bandname
if not opts.skip_rename_master:
    dstr = os.path.basename(fnams[0])[0:8]
    try:
        dtim = datetime.strptime(dstr,'%Y%m%d')
    except Exception:
        raise ValueError('Error in filename >>> '+fnams[0])
    band_list = list(data_1.getBandNames())
    bands = []
    for band in [band_list[j] for j in opts.band]:
        if dstr in band:
            continue
        band_new = band+'_'+dstr
        ProductUtils.copyBand(band,data_1,band_new,data_1,True)
        bands.append(band_new)
    if len(bands) > 0:
        params = HashMap()
        params.put('sourceBands',','.join(bands))
        data_tmp = GPF.createProduct('BandSelect',params,data_1)
        data_1 = data_tmp
# Collocation
for i in range(1,len(fnams)):
    # Read original product
    data_2 = ProductIO.readProduct(fnams[i])
    # Attach bandname
    if not opts.skip_rename_slave:
        dstr = os.path.basename(fnams[i])[0:8]
        try:
            dtim = datetime.strptime(dstr,'%Y%m%d')
        except Exception:
            raise ValueError('Error in filename >>> '+fnams[i])
        band_list = list(data_2.getBandNames())
        bands = []
        for band in [band_list[j] for j in opts.band]:
            if dstr in band:
                continue
            band_new = band+'_'+dstr
            ProductUtils.copyBand(band,data_2,band_new,data_2,True)
            bands.append(band_new)
        if len(bands) > 0:
            params = HashMap()
            params.put('sourceBands',','.join(bands))
            data_tmp = GPF.createProduct('BandSelect',params,data_2)
            data_2 = data_tmp
    # Collocate data
    params = HashMap()
    params.put('ResamplingType','Nearest neighbour')
    params.put('renameMasterComponents',False)
    params.put('renameSlaveComponents',False)
    products = HashMap()
    products.put("master",data_1)
    products.put("slave",data_2)
    data_tmp = GPF.createProduct('Collocate',params,products)
    data_1 = data_tmp
# Remove flag bands
band_list = list(data_1.getBandNames())
bands = []
for band in band_list:
    if not re.search('flag',band.lower()):
        bands.append(band)
if len(bands) != len(band_list):
    params = HashMap()
    params.put('sourceBands',','.join(bands))
    data_tmp = GPF.createProduct('BandSelect',params,data_1)
    data_1 = data_tmp
# Output image
ProductIO.writeProduct(data_1,opts.output_fnam,'GeoTiff')
