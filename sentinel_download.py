#!/usr/bin/env python
import os
import psutil
import sys
import re
import time
import signal
import atexit
import numpy as np
from datetime import datetime
from glob import glob
from subprocess import Popen,check_output,PIPE
from sentinelsat import SentinelAPI
from optparse import OptionParser,IndentedHelpFormatter

TIMEOUT = 60
DOWNLOAD_CHECK_TIME = 10
ONLINE_CHECK_TIME = 300
WAIT_TIME = 10
MAX_RETRY = 100

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
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
parser.add_option('-P','--path',default=None,help='Set the path where the files will be saved.')
parser.add_option('-q','--query',default=None,help='Extra search keywords you want to use in the query. Separate keywords with comma. Example: \'producttype=GRD,polarisationmode=HH\'.')
parser.add_option('-T','--timeout',default=TIMEOUT,type='float',help='Timeout to download data in sec (%default)')
parser.add_option('-w','--download_check_time',default=DOWNLOAD_CHECK_TIME,type='int',help='Wait time to check downloaded data in sec (%default)')
parser.add_option('-W','--wait_time',default=WAIT_TIME,type='int',help='Wait time to download data in sec (%default)')
parser.add_option('-O','--online_check_time',default=ONLINE_CHECK_TIME,type='int',help='Wait time to check online data in sec (%default)')
parser.add_option('-M','--max_retry',default=MAX_RETRY,type='int',help='Maximum number of retries to download data (%default)')
parser.add_option('-d','--download',default=False,action='store_true',help='Download all results of the query. (%default)')
parser.add_option('-C','--checksum',default=False,action='store_true',help='Verify the downloaded files\' integrity by checking its MD5 checksum. (%default)')
parser.add_option('-f','--footprints',default=False,action='store_true',help='Create a geojson file search_footprints.geojson with footprints and metadata of the returned products. (%default)')
parser.add_option('-v','--version',default=False,action='store_true',help='Show the version and exit. (%default)')
parser.add_option('-V','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
(opts,args) = parser.parse_args()

process_id = None
def kill_process():
    if process_id is not None:
        if not psutil.pid_exists(process_id):
            return
        sys.stderr.write('KILL PROCESS\n')
        if os.name.lower() != 'posix':
            parent = psutil.Process(process_id)
            children = parent.children(recursive=True)
            for child in children:
                child.kill()
            parent.kill()
        else:
            os.killpg(os.getpgid(process_id),signal.SIGTERM)
atexit.register(kill_process)

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
#if opts.download:
#    command += ' --download'
if not opts.checksum:
    command += ' --no-checksum'
if opts.footprints:
    command += ' --footprints'
if opts.version:
    command += ' --version'
sys.stderr.write(command+'\n')
out = check_output(command+' 2>&1',shell=True).decode()
sys.stderr.write(out+'\n')
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
if opts.url is not None:
    api = SentinelAPI(opts.user,opts.password,opts.url)
else:
    api = SentinelAPI(opts.user,opts.password)
for i,uuid in enumerate(uuids):
    out = api.get_product_odata(uuid)
    name = out['title']
    size = out['size']
    stat = out['Online']
    names.append(name)
    sizes.append(size)
    stats.append(stat)
    sys.stderr.write('{:4d} {:40s} {:70s} {:10d} {:7s}\n'.format(i+1,uuid,name,size,'Online' if stat else 'Offline'))
api.session.close() # has any effect?

if opts.download:
    path = '.' if opts.path is None else opts.path
    for i in range(len(uuids)):
        # Check data availability
        fnam = os.path.join(path,names[i]+'.zip')
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
            time.sleep(opts.online_check_time)
            continue
        api.session.close() # has any effect?
        # Download data
        gnam = os.path.join(fnam+'.incomplete')
        command = 'sentinelsat'
        command += ' --download'
        command += ' --uuid {}'.format(uuids[i])
        if opts.user is not None:
            command += ' --user {}'.format(opts.user)
        if opts.password is not None:
            command += ' --password {}'.format(opts.password)
        if opts.url is not None:
            command += ' --url {}'.format(opts.url)
        if opts.path is not None:
            command += ' --path {}'.format(opts.path)
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
                        break
                    else:
                        os.rename(fnam,gnam)
            # Start a process
            p = None
            process_id = None
            start_time = time.time()
            try:
                if opts.verbose:
                    if os.name.lower() != 'posix':
                        p = Popen(command,shell=True)
                    else:
                        p = Popen(command,shell=True,preexec_fn=os.setsid)
                else:
                    if os.name.lower() != 'posix':
                        p = Popen(command,stdout=PIPE,stderr=PIPE,bufsize=1,universal_newlines=True,shell=True)
                    else:
                        p = Popen(command,stdout=PIPE,stderr=PIPE,bufsize=1,universal_newlines=True,shell=True,preexec_fn=os.setsid)
                process_id = p.pid
            except Exception:
                sys.stderr.write('Failed to run the command. Wait for {} sec\n'.format(opts.wait_time))
                time.sleep(opts.wait_time)
                continue
            while True: # loop to terminate the process
                if not opts.verbose:
                    for line in p.stderr.readline().splitlines():
                        line_lower = line.lower()
                        if re.search('raise ',line_lower):
                            continue
                        elif re.search('sentinelapiltaerror',line_lower):
                            continue
                        elif re.search('will ',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('triggering ',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('requests ',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('accept',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('quota',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('downloading ',line_lower):
                            sys.stderr.write(line+'\n')
                        elif re.search('downloading:',line_lower):
                            sys.stderr.write(line+'\r')
                    sys.stderr.flush()
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
            kill_process()
            if not opts.verbose:
                out,err = p.communicate()
                sys.stderr.write('\n')
                for line in err.decode().splitlines():
                    line_lower = line.lower()
                    if re.search('raise ',line_lower):
                        continue
                    elif re.search('sentinelapiltaerror',line_lower):
                        continue
                    elif re.search('will ',line_lower):
                        sys.stderr.write(line+'\n')
                    elif re.search('triggering ',line_lower):
                        sys.stderr.write(line+'\n')
                    elif re.search('requests ',line_lower):
                        sys.stderr.write(line+'\n')
                    elif re.search('accept',line_lower):
                        sys.stderr.write(line+'\n')
                    elif re.search('quota',line_lower):
                        sys.stderr.write(line+'\n')
                    elif re.search('downloading ',line_lower):
                        sys.stderr.write(line+'\n')
            process_id = None
            sys.stderr.write('Wait for {} sec\n'.format(opts.wait_time))
            time.sleep(opts.wait_time)
