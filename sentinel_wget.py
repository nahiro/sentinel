#!/usr/bin/env python
import os
import psutil
import sys
import re
import time
import numpy as np
import xml.etree.ElementTree as ET
from collections import OrderedDict
from datetime import datetime
from glob import glob
from subprocess import Popen,check_output,PIPE,STDOUT
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
SERVER = 'https://scihub.copernicus.eu/dhus'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('--server',default=SERVER,help='Server URL (%default)')
parser.add_option('-u','--user',default=None,help='Username (or environment variable DHUS_USER is set)')
parser.add_option('-p','--password',default=None,help='Password (or environment variable DHUS_PASSWORD is set)')
parser.add_option('-U','--url',default=None,help='API URL (or environment variable DHUS_URL is set) (%default)')
parser.add_option('-s','--start',default=None,help='Start date of the query in the format YYYYMMDD.')
parser.add_option('-e','--end',default=None,help='End date of the query in the format YYYYMMDD.')
parser.add_option('-g','--geometry',default=None,help='Search area geometry as GeoJSON file.')
parser.add_option('-I','--uuid',default=None,help='Select a specific product UUID instead of a query. Multiple UUIDs can be separated by comma.')
parser.add_option('-N','--name',default=None,help='Select specific product(s) by filename. Multiple names can be separated by comma. Supports wildcards.')
parser.add_option('-S','--sentinel',default=None,help='[1|2|3|5] Limit search to a Sentinel satellite (constellation)')
parser.add_option('-i','--instrument',default=None,help='[MSI|SAR-C SAR|SLSTR|OLCI|SRAL] Limit search to a specific instrument on a Sentinel satellite.')
parser.add_option('-t','--producttype',default=None,help='Limit search to a Sentinel product type.')
parser.add_option('-c','--cloud',default=None,type='int',help='Maximum cloud cover in percent. (requires --sentinel to be 2 or 3)')
parser.add_option('-o','--order_by',default=None,help='Comma-separated list of keywords to order the result by. Prefix keywords with \'-\' for descending order.')
parser.add_option('-l','--limit',default=None,type='int',help='Maximum number of results to return.  Defaults to no limit.')
parser.add_option('-L','--log',default=None,help='Set the log file name.')
parser.add_option('-P','--path',default=None,help='Set the path where the files will be saved.')
parser.add_option('-q','--query',default=None,help='Extra search keywords you want to use in the query. Separate keywords with comma. Example: \'producttype=GRD,polarisationmode=HH\'.')
parser.add_option('-d','--download',default=False,action='store_true',help='Download all results of the query. (%default)')
parser.add_option('-Y','--sort_year',default=False,action='store_true',help='Sort files by year. (%default)')
parser.add_option('-C','--checksum',default=False,action='store_true',help='Verify the downloaded files\' integrity by checking its MD5 checksum. (%default)')
parser.add_option('-f','--footprints',default=False,action='store_true',help='Create a geojson file search_footprints.geojson with footprints and metadata of the returned products. (%default)')
parser.add_option('-v','--version',default=False,action='store_true',help='Show the version and exit. (%default)')
parser.add_option('-Q','--quiet',default=False,action='store_true',help='Quiet mode (%default)')
(opts,args) = parser.parse_args()

command = 'sentinelsat'
if opts.user is not None:
    command += ' --user {}'.format(opts.user)
if opts.password is not None:
    command += ' --password {}'.format(opts.password)
if opts.url is not None:
    command += ' --url {}'.format(opts.url)
if opts.start is not None:
    command += ' --start {}'.format(opts.start)
if opts.end is not None:
    command += ' --end {}'.format(opts.end)
if opts.geometry is not None:
    command += ' --geometry {}'.format(opts.geometry)
if opts.uuid is not None:
    command += ' --uuid {}'.format(opts.uuid)
if opts.name is not None:
    command += ' --name {}'.format(opts.name)
if opts.sentinel is not None:
    command += ' --sentinel {}'.format(opts.sentinel)
if opts.instrument is not None:
    command += ' --instrument {}'.format(opts.instrument)
if opts.producttype is not None:
    command += ' --producttype {}'.format(opts.producttype)
if opts.cloud is not None:
    command += ' --cloud {}'.format(opts.cloud)
if opts.order_by is not None:
    command += ' --order_by {}'.format(opts.order_by)
if opts.limit is not None:
    command += ' --limit {}'.format(opts.limit)
if opts.path is not None:
    command += ' --path {}'.format(opts.path)
if opts.query is not None:
    command += ' --query {}'.format(opts.query)
if not opts.checksum:
    command += ' --no-checksum'
if opts.footprints:
    command += ' --footprints'
if opts.version:
    command += ' --version'
sys.stderr.write(command+'\n')
sys.stderr.flush()
#out = check_output(command,shell=True,stderr=PIPE).decode()
out = check_output(command,shell=True,stderr=STDOUT).decode()
sys.stderr.write(out+'\n')
sys.stderr.flush()
if opts.version:
    sys.exit()

uuids = []
for line in out.splitlines():
    # Product bb7a7783-f91b-4a35-907d-a6ddb807da73 - Date: 2019-07-09T02:55:59.024Z, Instrument: MSI, Mode: , Satellite: Sentinel-2, Size: 1.09 GB
    m = re.search('Product\s+(\S+)',line)
    if not m:
        continue
    uuids.append(m.group(1))

names = []
sizes = []
stats = []
for i,uuid in enumerate(uuids):
    command = 'wget'
    command += ' --no-check-certificate'
    command += ' --output-document -'
    command += ' "'+opts.server+'/odata/v1/Products(\'{}\')"'.format(uuid)
    out = check_output(command,shell=True,stderr=PIPE).decode()
    root = ET.fromstring(out)
    child = None
    for value in root:
        prefix,has_namespace,postfix = value.tag.rpartition('}')
        if postfix.lower() == 'properties':
            child = value
            break
    tags = OrderedDict()
    for value in child:
        prefix,has_namespace,postfix = value.tag.rpartition('}')
        if len(value) > 0:
            tags[postfix] = OrderedDict()
            for value_2 in value:
                prefix_2,has_namespace_2,postfix_2 = value_2.tag.rpartition('}')
                tags[postfix][postfix_2] = value_2.text
        else:
            tags[postfix] = value.text
    name = tags['Name']
    size = int(tags['ContentLength'])
    stat = eval(tags['Online'].capitalize())
    names.append(name)
    sizes.append(size)
    stats.append(stat)
    sys.stderr.write('{:4d} {:40s} {:70s} {:10d} {:7s}\n'.format(i+1,uuid,name,size,'Online' if stat else 'Offline'))
    sys.stderr.flush()

"""
if opts.download:
    path = '.' if opts.path is None else opts.path
    for i in range(len(uuids)):
        # Check data availability
        if opts.sort_year:
            m = re.search('^\S+[^\d]_('+'\d'*4+')'+'\d'*4+'T'+'\d'*6+'_',names[i])
            if not m:
                raise ValueError('Error in file name >>> '+names[i])
            year = m.group(1)
            dnam = os.path.join(path,year)
            if not os.path.exists(dnam):
                os.makedirs(dnam)
            if not os.path.isdir(dnam):
                raise IOError('Error, no such directory >>> '+dnam)
        else:
            dnam = path
        fnam = os.path.join(dnam,names[i]+'.zip')
        gnam = os.path.join(fnam+'.incomplete')
        # Skip if gnam does not exist and fnam with expected size exists
        if os.path.exists(fnam) and not os.path.exists(gnam):
            fsiz = os.path.getsize(fnam)
            if fsiz == sizes[i]:
                sys.stderr.write('###### Successfully downloaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
                continue
        gnam = os.path.join(fnam+'.request')
        if opts.url is not None:
            api = SentinelAPI(opts.user,opts.password,opts.url)
        else:
            api = SentinelAPI(opts.user,opts.password)
        while True:
            out = api.get_product_odata(uuids[i])
            stat = out['Online']
            if stat: # Online
                if os.path.exists(gnam):
                    os.remove(gnam)
                break
            sys.stderr.write('Offline. Wait for {} sec >>> {}\n'.format(opts.online_check_time,fnam))
            sys.stderr.flush()
            time.sleep(opts.online_check_time)
            continue
        api.session.close() # has any effect?
        # Download data
        gnam = os.path.join(fnam+'.incomplete')
        command = 'sentinel_wget.py'
        command += ' --download'
        command += ' --uuid {}'.format(uuids[i])
        if opts.user is not None:
            command += ' --user {}'.format(opts.user)
        if opts.password is not None:
            command += ' --password {}'.format(opts.password)
        if opts.url is not None:
            command += ' --url {}'.format(opts.url)
        command += ' --path {}'.format(dnam)
        if not opts.checksum:
            command += ' --no-checksum'
        for ntry in range(opts.max_retry): # loop to download 1 file
            # Exit if gnam does not exist and fnam with expected size exists
            if os.path.exists(fnam):
                fsiz = os.path.getsize(fnam)
                if os.path.exists(gnam):
                    gsiz = os.path.getsize(gnam)
                    if gsiz > fsiz:
                        os.remove(fnam)
                    else:
                        os.remove(gnam)
                        os.rename(fnam,gnam)
                else:
                    if fsiz == sizes[i]:
                        sys.stderr.write('###### Successfully downloaded >>> {}\n'.format(fnam))
                        sys.stderr.flush()
                        if opts.log is not None:
                            with open(opts.log,'a') as fp:
                                fp.write(fnam+'\n')
                        break
                    else:
                        os.rename(fnam,gnam)
            # Start a process
            p = None
            process_id = None
            start_time = time.time()
            try:
                if opts.quiet:
                    if os.name.lower() != 'posix':
                        p = Popen(command,stdout=PIPE,stderr=PIPE,bufsize=0,universal_newlines=True,encoding='utf-8',shell=True)
                    else:
                        p = Popen(command,stdout=PIPE,stderr=PIPE,bufsize=0,universal_newlines=True,encoding='utf-8',shell=True,preexec_fn=os.setsid)
                else:
                    if os.name.lower() != 'posix':
                        p = Popen(command,universal_newlines=True,encoding='utf-8',shell=True)
                    else:
                        p = Popen(command,universal_newlines=True,encoding='utf-8',shell=True,preexec_fn=os.setsid)
                process_id = p.pid
            except Exception:
                sys.stderr.write('Failed to run the command. Wait for {} sec\n'.format(opts.wait_time))
                sys.stderr.flush()
                time.sleep(opts.wait_time)
                continue
            while True: # loop to terminate the process
                # Exit if fnam exists or gnam exists and its size does not change for opts.timeout seconds
                if p.poll() is not None: # the process has terminated
                    break
                elif os.path.exists(fnam):
                    break
                elif os.path.exists(gnam):
                    diff = time.time()-max(start_time,os.path.getmtime(gnam))
                    if diff > opts.timeout:
                        break
                time.sleep(opts.download_check_time)
            if p.poll() is None: # the process hasn't terminated yet.
                sys.stderr.write('\nTimeout\n')
                sys.stderr.flush()
            kill_process()
            if opts.quiet:
                out,err = p.communicate()
                sys.stderr.write('\n')
                sys.stderr.flush()
"""
