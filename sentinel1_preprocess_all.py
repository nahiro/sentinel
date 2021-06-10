#!/usr/bin/env python
import os
import sys
import re
from datetime import datetime,timedelta
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
SCRDIR = os.path.join(HOME,'Script')
DATDIR = os.path.join(HOME,'Work','Sentinel-1')
SITES = ['Bojongsoang','Cihea']

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=None,help='End date in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES
if opts.end is None:
    opts.end = '{:%Y%m%d}'.format(datetime.now())
dmax = datetime.strptime(opts.end,'%Y%m%d')
if opts.str is None:
    opts.str = '{:%Y%m%d}'.format(dmax-timedelta(days=1))
dmin = datetime.strptime(opts.str,'%Y%m%d')

for site in opts.sites:
    datdir = os.path.join(opts.datdir,site)
    fnams = []
    dstrs = []
    for d in sorted(os.listdir(datdir)):
        if not re.search('^\d\d\d\d$',d):
            continue
        year = int(d)
        if year < dmin.year or year > dmax.year:
            continue
        dnam = os.path.join(datdir,d)
        for f in sorted(os.listdir(dnam)):
            #S1A_IW_GRDH_1SDV_20171227T223338_20171227T223405_019894_021DC8_434F.zip
            #S1B_IW_GRDH_1SDV_20200116T223300_20200116T223336_019848_025883_2DEF.zip
            m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*8+')T\S+\.zip$',f)
            if not m:
                continue
            dstr = m.group(1)
            dtim = datetime.strptime(dstr,'%Y%m%d')
            if (dtim < dmin) or (dtim > dmax):
                continue
            fnam = os.path.join(dnam,f)
            fnams.append(fnam)
            dstrs.append(dstr)
    if len(dstrs) < 1:
        continue

    for fnam in fnams:
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'sentinel1_preprocess.py')
        command += ' '+fnam
        command += ' --datdir '+os.path.join(datdir,'sigma0_speckle')
        command += ' --site '+site
        command += ' --speckle'
        command += ' --iangle_value'
        command += ' --std_grid'
        command += ' --tiff'
        call(command,shell=True)
        # Remove cache
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'remove_snap_cache.py')
        call(command,shell=True)
    for dstr in dstrs:
        fnam = os.path.join(datdir,'sigma0_speckle',dstr+'.tif')
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'sentinel_resample.py')
        command += ' '+fnam
        command += ' --datdir '+os.path.join(datdir,'sigma0_speckle')
        command += ' --site '+site
        command += ' --read_comments'
        call(command,shell=True)
