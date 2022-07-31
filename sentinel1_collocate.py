# SENTINEL-1 data Preprocessing -------------------------------------------------------#
import os
import sys
import glob
os.environ['_JAVA_OPTIONS'] = '-Xmx16240m'     # Seve memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx16240m')   # Seve memory for JAVA
import snappy
from snappy import jpy
from snappy import ProductIO, WKTReader
from snappy import GPF
from snappy import Product
from snappy import ProductData
from snappy import ProductUtils
from snappy import HashMap
import numpy as np
try:
    import gdal
except Exception:
    from osgeo import gdal

# Collocation ------------------------------------------------------------------------#
def collocate(input_img1, input_img2, write_product=False, write_fnam=None):
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    # Read data
    if type(input_img1) is str:
        data_master = ProductIO.readProduct(input_img1)
    else:
        data_master = input_img1
    if type(input_img2) is str:
        data_slave = ProductIO.readProduct(input_img2)
    else:
        data_master = input_img2
    parameters = HashMap()
    parameters.put('ResamplingMethod', 'Nearest neighbour')
    sourceProducts = HashMap()
    sourceProducts.put("master", data_master)
    sourceProducts.put("slave", data_slave)
    output_img = GPF.createProduct('Collocate', parameters, sourceProducts)
    if write_product:
        if write_fnam is None:
            fnam = 'collocate.tif'
        else:
            fnam = write_fnam
        ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
    return output_img

def collocate_all(input_imgs, write_product):
    if len(input_imgs) < 2:
        raise ValueError('Error, len(input_imgs)={}'.format(len(input_imgs)))
    col = collocate(input_imgs[0],input_imgs[1])
    for i in range(2,len(input_imgs)):
        col = collocate(col,input_imgs[i])
    if write_product:
        collocation = "collocation_all.tif"
        ProductIO.writeProduct(col, collocation, 'GeoTiff')
    return col

datdir = '../2018'
write_flag = True
polarization = 'VH'

# Create a filelist
filelist_safe = glob.glob(os.path.join(datdir,'*.SAFE'))

data = []  # Empty list
for f in sorted(filelist_safe):
    date = os.path.basename(f).split('_')[4][0:8]
    fnam = date+'_'+polarization+'_dB.dim'
    if os.path.exists(fnam):
        data.append(fnam)
collocate_all(data,write_flag)
