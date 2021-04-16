#!/usr/bin/env python
import os
import sys
import shutil
import re
from datetime import datetime
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
TOPDIR = '/mnt/hlab/Data/Transplanting_date'
SITE = 'Cihea'
LEVEL = 'test'
VERSION = 'v1.0'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--topdir',default=TOPDIR,help='Top directory (%default)')
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('--site',default=SITE,help='Site name (%default)')
parser.add_option('--level',default=LEVEL,help='Analysis level, test/final/preliminary (%default)')
parser.add_option('--version',default=VERSION,help='Product version (%default)')
parser.add_option('--date',default=None,help='Date in the format YYYYMMDD (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
if opts.date is None:
    opts.date = datetime.now().strftime('%Y%m%d')
fnams = args

# Destination directory
dstr_year = datetime.strptime(opts.date,'%Y%m%d').strftime('%Y')
dstdir = os.path.join(opts.topdir,opts.site,opts.level,opts.version,dstr_year,opts.date)
if not os.path.exists(dstdir):
    os.makedirs(dstdir)
if not os.path.isdir(dstdir):
    raise IOError('Error, no such directory >>> '+dstdir)

for input_fnam in fnams:
    copy_fnam = os.path.basename(input_fnam)
    flag = True
    for f in sorted(os.listdir(dstdir)):
        if f.upper() == copy_fnam.upper():
            if opts.overwrite:
                sys.stderr.write('File exists, delete >>> '+f+'\n')
                os.remove(os.path.join(dstdir,f))
            else:
                sys.stderr.write('File exists, skip   >>> '+f+'\n')
                flag = False # no need to upload
                break
    if flag: # copy file
        shutil.copy2(input_fnam,os.path.join(dstdir,copy_fnam))
        sys.stderr.write('Successfully copied >>> '+copy_fnam+'\n')
