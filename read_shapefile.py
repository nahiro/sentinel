#!/usr/bin/env python
import sys
import shapefile
from optparse import OptionParser,IndentedHelpFormatter

# Default values
LENG = 20

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-l','--leng',default=LENG,type='int',help='Item length (%default)')
parser.add_option('-n','--field_name',default=False,action='store_true',help='Output field name only (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
inp_fnam = args[0]

item_format = ' {{:>{:d}s}}'.format(opts.leng)

r = shapefile.Reader(inp_fnam)
sys.stdout.write('{:>6s}'.format('Index'))
for rec in r.fields[-len(r.record()):]:
    sys.stdout.write(item_format.format(str(rec[0]).strip()))
sys.stdout.write('\n')
if opts.field_name:
    sys.exit()
for i,shaperec in enumerate(r.iterShapeRecords()):
    sys.stdout.write('{:6d}'.format(i))
    for rec in shaperec.record:
        sys.stdout.write(item_format.format(str(rec).strip()))
    sys.stdout.write('\n')

