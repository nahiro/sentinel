#!/usr/bin/env python
import os
import sys
import shutil
import re
from datetime import datetime,timedelta
import numpy as np
import shapefile
import cartopy.io.shapereader as shpreader
from matplotlib.dates import date2num,num2date
from optparse import OptionParser,IndentedHelpFormatter

# Default values
FIELD_NAME = 'trans_t1'
FIELD_TYPE = 'character'
FIELD_LENGTH = 13
DECIMAL_LENGTH = 0
OUT_FNAM = 'output'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-a','--field_name',default=FIELD_NAME,help='Field name to add (%default)')
parser.add_option('-t','--field_type',default=FIELD_TYPE,help='Field type, character, float, or number (%default)')
parser.add_option('-l','--field_length',default=FIELD_LENGTH,type='int',help='Field length (%default)')
parser.add_option('-d','--decimal_length',default=DECIMAL_LENGTH,type='int',help='Decimal length (%default)')
parser.add_option('-s','--shp_fnam',default=None,help='Input shapefile name (%default)')
parser.add_option('-o','--out_fnam',default=OUT_FNAM,help='Output shapefile name (%default)')
parser.add_option('-O','--output_field_name',default=None,action='append',help='Output field name (%default)')
parser.add_option('-T','--num2date',default=False,action='store_true',help='Convert number to text date (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args
field_type = opts.field_type[0].upper()
if not field_type.upper() in ['C','F','N']:
    raise ValueError('Error in field type >>> '+opts.field_type)

r = shapefile.Reader(opts.shp_fnam)
nobject = len(r)

data = []
if opts.output_field_name is None:
    dtim = []
for fnam in fnams:
    dtmp = []
    records = list(shpreader.Reader(fnam).records())
    for rec in records:
        dtmp.append(rec.attributes[opts.field_name])
    if len(dtmp) != nobject:
        raise ValueError('Error, len(dtmp)={}, nobject={}'.format(len(dtmp),nobject))
    if opts.output_field_name is None:
        f = os.path.basename(fnam)
        m = re.search('\D('+'\d'*8+')\D',f)
        if not m:
            m = re.search('^('+'\d'*8+')\D',f)
            if not m:
                raise ValueError('Error in finding date >>> '+f)
        dstr = m.group(1)
        d = datetime.strptime(dstr,'%Y%m%d')
        dtim.append(d)
    data.append(dtmp)
data = np.array(data)
ndat = len(data)
if opts.output_field_name is None:
    dtim = np.array(dtim)
    if len(dtim) != ndat:
        raise ValueError('Error, len(dtim)={}, ndat={}'.format(len(dtim),ndat))

w = shapefile.Writer(opts.out_fnam)
w.shapeType = shapefile.POLYGON
w.fields = r.fields[1:] # skip first deletion field
if opts.output_field_name is None:
    for i in range(ndat):
        w.field('{:%Y/%m/%d}'.format(dtim[i]),field_type,opts.field_length,opts.decimal_length)
else:
    if len(opts.output_field_name) != ndat:
        raise ValueError('Error, len(opts.output_field_name)={}, ndat={}'.format(len(opts.output_field_name),ndat))
    for i in range(ndat):
        w.field(opts.output_field_name[i],field_type,opts.field_length,opts.decimal_length)
for iobj,shaperec in enumerate(r.iterShapeRecords()):
    rec = shaperec.record
    shp = shaperec.shape
    if field_type == 'F':
        data_list = [float(d) for d in data[:,iobj]]
    elif field_type == 'N':
        data_list = [int(d) for d in data[:,iobj]]
    else:
        if opts.num2date:
            data_list = ['N/A' if np.isnan(d) else num2date(np.round(d)+0.1).strftime('%Y/%m/%d') for d in data[:,iobj]]
        else:
            data_list = list(data[:,iobj])
    rec.extend(data_list)
    w.shape(shp)
    w.record(*rec)
w.close()
shutil.copy2(opts.shp_fnam+'.prj',opts.out_fnam+'.prj')
