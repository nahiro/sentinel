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
L2ADIR = os.path.join(HOME,'Work','Sentinel-2','L2A')
SITES = ['Bojongsoang','Cihea']
SUBDIRS = ['Cihea:sigma0','Bojongsoang:sigma0_speckle']
DATE_FINAL = 5
MAX_RETRY = 10

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--scrdir',default=SCRDIR,help='Script directory (%default)')
parser.add_option('--datdir',default=DATDIR,help='Data directory (%default)')
parser.add_option('--l2adir',default=L2ADIR,help='L2A data directory (%default)')
parser.add_option('-s','--str',default=None,help='Start date in the format YYYYMMDD (%default)')
parser.add_option('-e','--end',default=None,help='End date in the format YYYYMMDD (%default)')
parser.add_option('-S','--sites',default=None,action='append',help='Target sites ({})'.format(SITES))
parser.add_option('--subdirs',default=None,action='append',help='Sub data directory, for example, Cihea:sigma0 ({})'.format(SUBDIRS))
parser.add_option('--date_final',default=DATE_FINAL,type='int',help='Date to calculate final estimation (%default)')
parser.add_option('-M','--max_retry',default=MAX_RETRY,type='int',help='Maximum number of retries to download data (%default)')
parser.add_option('--skip_upload',default=False,action='store_true',help='Skip upload (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
(opts,args) = parser.parse_args()
if opts.sites is None:
    opts.sites = SITES
if opts.subdirs is None:
    opts.subdirs = SUBDIRS
subdir = {}
for s in opts.subdirs:
    m = re.search('([^:]+):([^:]+)',s)
    if not m:
        raise ValueError('Error in subdir >>> '+s)
    subdir.update({m.group(1).lower():m.group(2)})

for site in opts.sites:
    site_low = site.lower()
    datdir = os.path.join(opts.datdir,site)
    log = os.path.join(datdir,site.lower()+'.log')
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'sentinel1_update.py')
    command += ' --scrdir '+opts.scrdir
    command += ' --datdir '+opts.datdir
    if opts.str is not None:
        command += ' --str '+opts.str
    if opts.end is not None:
        command += ' --end '+opts.end
    if opts.skip_upload:
        command += ' --skip_upload'
    command += ' --sites '+site
    #sys.stderr.write(command+'\n')
    call(command,shell=True)
    fnams = []
    dstrs = []
    if os.path.exists(log):
        with open(log,'r') as fp:
            for line in fp:
                item = line.split()
                if len(item) < 1:
                    continue
                fnam = item[0]
                f = os.path.basename(fnam)
                #S1A_IW_GRDH_1SDV_20171227T223338_20171227T223405_019894_021DC8_434F.zip
                #S1B_IW_GRDH_1SDV_20200116T223300_20200116T223336_019848_025883_2DEF.zip
                m = re.search('^S1[AB]_IW_GRDH_1SDV_('+'\d'*8+')T\S+\.zip$',f)
                if not m:
                    continue
                fnams.append(fnam)
                dstrs.append(m.group(1))
    if len(dstrs) < 1:
        continue
    for fnam in fnams:
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'sentinel1_preprocess.py')
        command += ' '+fnam
        command += ' --site '+site
        command += ' --datdir '+os.path.join(datdir,subdir[site_low])
        if 'speckle' in subdir[site_low].lower():
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
        fnam = os.path.join(datdir,subdir[site_low],dstr+'.tif')
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'sentinel_resample.py')
        command += ' '+fnam
        command += ' --datdir '+os.path.join(datdir,subdir[site_low])
        command += ' --site '+site
        command += ' --read_comments'
        call(command,shell=True)
        dtim = datetime.strptime(dstr,'%Y%m%d')
        d1 = (dtim+timedelta(days=-1)).strftime('%Y%m%d')
        d2 = (dtim+timedelta(days=+1)).strftime('%Y%m%d')
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'get_preliminary_estimation.py')
        command += ' --scrdir '+opts.scrdir
        command += ' --datdir '+opts.datdir
        command += ' --str '+d1
        command += ' --end '+d2
        command += ' --sites '+site
        call(command,shell=True)
    if os.path.exists(log):
        os.remove(log)

dcur = datetime.now()
if dcur.day == opts.date_final:
    for site in opts.sites:
        command = 'python'
        command += ' '+os.path.join(opts.scrdir,'get_final_estimation.py')
        command += ' --scrdir '+opts.scrdir
        command += ' --datdir '+opts.datdir
        command += ' --sites '+site
        call(command,shell=True)

for site in ['Bojongsoang']:
    datdir = os.path.join(opts.datdir,site)
    command = 'python'
    command += ' '+os.path.join(opts.scrdir,'sentinel2_update.py')
    command += ' --scrdir '+opts.scrdir
    command += ' --datdir '+opts.l2adir
    if opts.str is not None:
        command += ' --str '+opts.str
    else:
        command += ' --str {:%Y%m%d}'.format(dcur+timedelta(days=-60))
    if opts.end is not None:
        command += ' --end '+opts.end
    else:
        command += ' --end {:%Y%m%d}'.format(dcur)
    if opts.skip_upload:
        command += ' --skip_upload'
    command += ' --sites '+site
    command += ' --max_retry {}'.format(opts.max_retry)
    #sys.stderr.write(command+'\n')
    call(command,shell=True)
