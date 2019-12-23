#!/usr/bin/env python
import os
import psutil
mem_size = int(psutil.virtual_memory().available*0.8e-6)
os.environ['_JAVA_OPTIONS'] = '-Xmx{}m'.format(mem_size)     # Save memory for JAVA
os.system('export _JAVA_OPTIONS=-Xmx{}m'.format(mem_size))   # Save memory for JAVA
import sys
from snappy import ProductIO
from optparse import OptionParser,IndentedHelpFormatter

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog filename [options]')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args[:]

for fnam in fnams:
    # Read original product
    data = ProductIO.readProduct(os.path.abspath(fnam))
    band_list = list(data.getBandNames())
    for i,band in enumerate(band_list):
        sys.stdout.write(os.path.basename(fnam)+' '+'{:4d}'.format(i+1)+' '+band+'\n')
