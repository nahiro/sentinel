#!/usr/bin/env python
import os
import sys
import shutil
import re
from datetime import datetime
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
TOPDIR = '/mnt/hlab/Data/Sentinel-2/L2A'
SITE = 'Bojongsoang'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('--topdir',default=TOPDIR,help='Top directory (%default)')
parser.add_option('-S','--site',default=SITE,help='Site name (%default)')
parser.add_option('-O','--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
parser.add_option('-D','--dry_run',default=False,action='store_true',help='Test mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args

for input_fnam in fnams:
    # S2A_MSIL2A_20210104T030121_N0214_R032_T48MYT_20210104T062157.zip
    fnam = os.path.basename(input_fnam)
    unam = fnam.upper()
    bnam,enam = os.path.splitext(unam)
    if enam != '.ZIP':
        if enam != '.SAFE':
            continue
    m = re.search('([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)',bnam)
    if not m:
        continue
    if m.group(1) != 'S2A' and m.group(1) != 'S2B':
        continue
    if m.group(2) != 'MSIL2A':
        sys.stderr.write('Warning, skipping file >>> '+fnam+'\n')
        sys.stderr.flush()
        continue
    d1 = datetime.strptime(m.group(3),'%Y%m%dT%H%M%S')
    d2 = datetime.strptime(m.group(7),'%Y%m%dT%H%M%S')
    if d1.date() != d2.date():
        sys.stderr.write('Warning, d1={}, d2={} >>> {}\n'.format(m.group(3),m.group(7),fnam))
        sys.stderr.flush()
    search_key = m.group(1)+'_'+m.group(2)+'_'+d1.strftime('%Y%m%d')
    copy_fnam = bnam+enam.lower()
    dstr_year = d1.strftime('%Y')
    dstdir = os.path.join(opts.topdir,opts.site,dstr_year) # Destination directory
    if not os.path.exists(dstdir):
        os.makedirs(dstdir)
    if not os.path.isdir(dstdir):
        raise IOError('Error, no such directory >>> '+dstdir)
    flag = True
    for f in sorted(os.listdir(dstdir)):
        if not re.search(search_key,f.upper()):
            continue
        if f.upper() == unam:
            if opts.overwrite:
                sys.stderr.write('File exists, delete   >>> '+f+'\n')
                sys.stderr.flush()
                os.remove(os.path.join(dstdir,f))
            else:
                sys.stderr.write('File exists, skip     >>> '+f+'\n')
                sys.stderr.flush()
                flag = False # no need to copy
                break
        else:
            sys.stderr.write('Warning, different file for the same date >>> '+f+'\n')
            sys.stderr.flush()
    if flag: # copy file
        if opts.dry_run:
            sys.stderr.write('cp {} {}\n'.format(input_fnam,os.path.join(dstdir,copy_fnam)))
            sys.stderr.flush()
        else:
            if os.path.isdir(input_fnam):
                shutil.copytree(input_fnam,os.path.join(dstdir,copy_fnam))
            else:
                shutil.copy2(input_fnam,os.path.join(dstdir,copy_fnam))
            sys.stderr.write('Successfully copied >>> '+copy_fnam+'\n')
            sys.stderr.flush()
