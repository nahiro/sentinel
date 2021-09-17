#!/usr/bin/env python
import sys
import numpy as np
import shapefile
from shapely.geometry import Point,Polygon
from optparse import OptionParser,IndentedHelpFormatter

# Default values
BLOCK_FNAM = '/home/naohiro/Work/SATREPS/Shapefile/cihea_testsite_200819/cihea_testsite_200819.shp'
FIELD_FNAM = '/home/naohiro/Work/SATREPS/Shapefile/field_GIS/cihea/New_Test_Sites.shp'
BLOCK_FIELD = 'Blok'
BLOCK = 'None'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-b','--block_fnam',default=BLOCK_FNAM,help='Block shape file name (%default)')
parser.add_option('-f','--field_fnam',default=FIELD_FNAM,help='Field shape file name (%default)')
parser.add_option('--block_field',default=BLOCK_FIELD,help='Field name of block in block_fnam (%default)')
parser.add_option('-B','--block',default=BLOCK,help='Default block name (%default)')
parser.add_option('--use_index',default=False,action='store_true',help='Use index instead of OBJECTID (%default)')
parser.add_option('--use_objectid',default=False,action='store_true',help='Use OBJECTID for default block name (%default)')
(opts,args) = parser.parse_args()

r1 = shapefile.Reader(opts.block_fnam)
blocks = []
polygons = []
for shaperec in r1.iterShapeRecords():
    rec = shaperec.record
    shp = shaperec.shape
    p1 = Polygon(shp.points)
    if not p1.is_valid:
        raise ValueError('Error, not valid >>> '+rec[1])
    blocks.append(getattr(rec,opts.block_field))
    polygons.append(p1)

r2 = shapefile.Reader(opts.field_fnam)
for iobj,shaperec in enumerate(r2.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    if opts.use_index:
        object_id = iobj+1
    else:
        object_id = getattr(rec,'OBJECTID')
    p2 = Polygon(shp.points)
    block = None
    ratios = []
    for i in range(len(blocks)):
        p1 = polygons[i]
        p3 = p2.intersection(p1)
        rat = p3.area/p2.area
        ratios.append(rat)
        if rat > 0.8:
            if block is None:
                block = blocks[i]
            else:
                raise ValueError('Error, block={}'.format(block))
    ratios = np.array(ratios)
    if block is None:
        i = np.argmax(ratios)
        if ratios[i] < 1.0e-10:
            if opts.use_objectid:
                block = str(object_id)
            else:
                block = opts.block
        else:
            block = blocks[i]
    sys.stdout.write('{:6d} {:>6s} {:.4e}\n'.format(object_id,block,ratios.max()))
    #break
