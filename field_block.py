#!/usr/bin/env python
import sys
import numpy as np
import shapefile
from shapely.geometry import Point,Polygon

shp_fnam1 = '/home/naohiro/Work/SATREPS/Shapefile/cihea_testsite_200819/cihea_testsite_200819.shp'
shp_fnam2 = '/home/naohiro/Work/SATREPS/Shapefile/field_GIS/cihea/New_Test_Sites.shp'

r1 = shapefile.Reader(shp_fnam1)
blocks = []
polygons = []
for shaperec in r1.iterShapeRecords():
    rec = shaperec.record
    shp = shaperec.shape
    blocks.append(rec[1])
    p1 = Polygon(shp.points)
    if not p1.is_valid:
        raise ValueError('Error, not valid >>> '+rec[1])
    polygons.append(p1)

r2 = shapefile.Reader(shp_fnam2)
for shaperec in r2.iterShapeRecords():
    rec = shaperec.record
    shp = shaperec.shape
    object_id = rec[0]
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
        block = blocks[i]
        if ratios[i] < 1.0e-10:
            block = 'None'
        sys.stderr.write('Warning, error in finding block for {} ({}, {})\n'.format(object_id,block,ratios[i]))
    sys.stdout.write('{:6d} {:>6s} {:.4e}\n'.format(object_id,block,1.0-ratios.max()))
