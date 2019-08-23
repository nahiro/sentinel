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
from osgeo import gdal

# Read original product
def readfile(fnam):
    sentinel = ProductIO.readProduct(fnam)
    band_list = list(sentinel.getBandNames())
    amp_VH = sentinel.getBand(band_list[0])
    w = amp_VH.getRasterWidth()
    h = amp_VH.getRasterHeight()
#   polarization = 'VH'
#   HashMap = snappy.jpy.get_type('java.util.HashMap')
#   jpy = snappy.jpy
    return sentinel,band_list,w,h


#Speckle-Filter
def speckle_filter(sentinel,polarization,write_product=False):
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()

    if type(polarization) is str:
        ps = [polarization]
    elif np.iterable(polarization):
        ps = polarization
    else:
        raise ValueError('Error in polarization >>> '+polarization)

    outputs = []
    for p in ps:
        if not p.upper in ['VH','VV']:
            raise ValueError('Error, unknown polarization >>>'+p)
        # Set the parameters
        parameters = HashMap()
        parameters.put('outputSigmaBand', True)
        parameters.put('sourceBands', 'Sigma0_' + p)
        parameters.put('Filter', 'Gamma Map')
        parameters.put('Filter Size X', 7)
        parameters.put('Filter Size Y', 7)
        output_img = GPF.createProduct('speckle_filter', parameters, sentinel)
        outputs.append(output_img)
        if write_product:
            if write_fnam is None:
                fnam = 'speckle_filter_'+p+'.dim'
            else:
                fnam = write_fnam+'_'+p+'.dim'
                ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
    if len(outputs) == 1:
        return outputs[0]
    else:
        return outputs



def calibration(sentinel,polarization,write_product=False, write_fnam=None):
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    filelist_safe = glob.glob('*.SAFE')
    if type(polarization) is str:  # One polarization or two polarizations
        ps = [polarization]
    elif np.iterable(polarization):
        ps = polarization
    else:
        raise ValueError('Error in polarization >>> '+polarization)
    outputs = []
    i = 0
    for p in ps:
        if not p.upper() in ['VH','VV']:   # Convert to an upper case
            raise ValueError('Error, unknown polarization >>>'+p)
        # Set the parameters
        parameters = HashMap()
        parameters.put('outputSigmaBand', True)
        parameters.put('sourceBands', 'Amplitude_' + p)
        parameters.put('selectedPolarisations', p)
        parameters.put('outputImageScaleInDb', False)
        output_img = GPF.createProduct('calibration', parameters, sentinel)
        outputs.append(output_img)
        if write_product:
            if write_fnam is None:
                fnam = 'calibration_'+p+'.dim'
            else:
                fnam = write_fnam+'_'+p+'.dim'
                ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
        i = i+1
    if len(outputs) == 1:
        return outputs[0]
    else:
        return outputs

# Subset image ------------------------------------------------------------------------#
def subset(input_img, write_product=False, write_fnam=None):
    if type(input_img) is str:
        data = ProductIO.readProduct(input_img)
    else:
        data = input_img
    SubsetOp = snappy.jpy.get_type('org.esa.snap.core.gpf.common.SubsetOp')
    WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader')
    wkt = "POLYGON((107.201 -6.898, 107.367 -6.898, 107.367 -6.767, 107.201 -6.767, 107.201 -6.898))"
    geom = WKTReader().read(wkt)
    # Set the parameters
    parameters = HashMap()
    parameters.put('copyMetadata', True)
    parameters.put('geoRegion', geom)
    parameters.put('outputImageScaleInDb', False)
    output_img = GPF.createProduct("Subset", parameters, data)
    if write_product:
        if write_fnam is None:
            fnam = 'subset.dim'
        else:
            fnam = write_fnam
        ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
    return output_img


# Terrain correction -----------------------------------------------------------------#
def terrain_correction(input_img, write_product=False, write_fnam=None):
    if type(input_img) is str:
        data = ProductIO.readProduct(input_img)
    else:
        data = input_img
    parameters = HashMap()
    parameters.put('demResamplingMethod', 'BILINEAR_INTERPOLATION')
    parameters.put('imgResamplingMethod', 'BILINEAR_INTERPOLATION')
    parameters.put('demName', 'SRTM 3Sec')
    parameters.put('pixelSpacingInMeter', 10.0)
    parameters.put('map projection', 'UTM/WGS84')
    parameters.put('sourceBands', 'Sigma0_VH')
    output_img = GPF.createProduct("Terrain-Correction", parameters, data)
    if write_product:
        if write_fnam is None:
            fnam = 'terrain_correction.tif'
        else:
            fnam = write_fnam
        ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
    return output_img

# Convert to dB values ---------------------------------------------------------------#
##dB = 10*log10(Amp)
def convert_dB(input_img, write_product=False, write_fnam=None):
    if type(input_img) is str:
        data = ProductIO.readProduct(input_img)
    else:
        data = input_img
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    BandDescriptor = jpy.get_type('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor')
    # Set the information of target band
    targetBand1 = BandDescriptor()
    targetBand1.name = 'Sigma0_VH'
    targetBand1.type = 'float32'
    targetBand1.expression = '10*log10(Sigma0_VH)'
    targetBands = jpy.array('org.esa.snap.core.gpf.common.BandMathsOp$BandDescriptor', 1)
    targetBands[0] = targetBand1
    parameters = HashMap()
    parameters.put('targetBands', targetBands)
    output_img = GPF.createProduct('BandMaths', parameters, data)
    if write_product:
        if write_fnam is None:
            fnam = 'convert_dB.tif'
        else:
            fnam = write_fnam
        ProductIO.writeProduct(output_img, fnam, 'BEAM-DIMAP')
    return output_img

datdir = '../../2018'
write_flag = False
polarization = 'VH'

# Create a filelist
filelist_safe = glob.glob(os.path.join(datdir,'*.SAFE'))

data = []  # Empty list
for f in sorted(filelist_safe):
    try:
        date = os.path.basename(f).split('_')[4][0:8]
        sentinel,band_list,w,h = readfile(f)
        print("f", f)
        print('subset')
        fnam = date+'_sub'
        target_0 = subset(sentinel, write_flag, fnam)
        fnam = date+'_calib'
        target_1 = calibration(target_0,polarization,write_flag,fnam)
        #print("target_1", target_1)
        print('terrain_correction')
        fnam = date+'_TC'
        target_2 = terrain_correction(target_1, write_flag, fnam)
        print('convert_dB')
        fnam = date+'_'+polarization+'_dB.dim'
        target_3 = convert_dB(target_2, write_product=True, write_fnam=fnam)
        #print("target_3", target_3)
    except Exception:
        continue
    data.append(fnam)
#print('collocate_all')
#filelist_dB = sorted(glob.glob(date+'_*_dB.dim'))
#collocate_all(data, True)
