#!/usr/bin/env python
import os
import sys
import re
import shutil
import zipfile
from glob import glob
from datetime import datetime,timedelta
import numpy as np
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

# Default values
DATDIR = '.'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-D','--datdir',default=DATDIR,help='Sentinel-1 data directory (%default)')
parser.add_option('-s','--start',default=None,help='Start date of the query in the format YYYYMMDD.')
parser.add_option('-e','--end',default=None,help='End date of the query in the format YYYYMMDD.')
parser.add_option('-g','--gamma0',default=False,action='store_true',help='Output gamma0 instead of sigma0 (%default)')
parser.add_option('--skip_orbit',default=False,action='store_true',help='Do not apply orbit file (%default)')
parser.add_option('--speckle',default=False,action='store_true',help='Apply speckle filter (%default)')
(opts,args) = parser.parse_args()
if opts.end is None:
    opts.end = datetime.now().strftime('%Y%m%d')
if opts.start is None:
    opts.start = (datetime.strptime(opts.end,'%Y%m%d')-timedelta(days=1)).strftime('%Y%m%d')

dmin = datetime.strptime(opts.start,'%Y%m%d')
dmax = datetime.strptime(opts.end,'%Y%m%d')

for year in range(dmin.year,dmax.year+1):
    flaglist = []
    filelist = []
    datalist = []
    datelist = []
    for f in glob(os.path.join(opts.datdir,str(year),'*.zip')):
        fnam = os.path.basename(f)
        m = re.search('^[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*.zip$',fnam)
        if m:
            dtim = datetime.strptime(m.group(1),'%Y%m%dT%H%M%S')
            if dtim < dmin or dtim > dmax:
                continue
            flaglist.append(False)
            filelist.append(f)
            datalist.append(f)
            datelist.append(dtim)
        else:
            g = None
            dtim = None
            z = zipfile.ZipFile(f)
            for d in z.namelist():
                dnam = os.path.dirname(d)
                if re.search('^\S+\.SAFE$',dnam):
                    m = re.search('^[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*.SAFE$',dnam)
                    if m:
                        g = os.path.splitext(dnam)[0]+'.zip'
                        dtim = datetime.strptime(m.group(1),'%Y%m%dT%H%M%S')
                        break
            if g is None:
                continue
            if not os.path.exists(g):
                flaglist.append(True)
            else:
                flaglist.append(False)
            filelist.append(f)
            datalist.append(os.path.abspath(g))
            datelist.append(dtim)
    flaglist = np.array(flaglist)
    filelist = np.array(filelist)
    datalist = np.array(datalist)
    datelist = np.array(datelist)
    indx = np.argsort(datelist)
    for flag,fnam,gnam,dtim in zip(flaglist[indx],filelist[indx],datalist[indx],datelist[indx]):
        dstr = dtim.strftime('%Y%m%d')
        sys.stderr.write(dstr+'\n')
        if os.path.exists('{}.dim'.format(dstr)):
            continue
        command = 'sentinel1_preprocess.py'
        command += ' '+gnam
        if opts.gamma0:
            command += ' --gamma0'
        if opts.skip_orbit:
            command += ' --skip_orbit'
        if opts.speckle:
            command += ' --speckle'
        if flag:
            os.symlink(fnam,gnam)
        call(command,shell=True)
        if flag:
            call('rm '+gnam,shell=True)
        call('rm -rf /tmp/snap-naohiro',shell=True)
        call('rm -rf /tmp/jffi*.tmp',shell=True)
        call('rm -rf /tmp/imageio*.tmp',shell=True)
        call('rm -rf /home/naohiro/.snap/var/cache/temp/imageio*.tmp',shell=True)
        #break
    #break
