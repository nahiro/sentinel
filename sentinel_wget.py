#!/usr/bin/env python
import os
import psutil
import sys
import re
import time
import hashlib
import xml.etree.ElementTree as ET
from collections import OrderedDict
from datetime import datetime
from subprocess import call,check_output,PIPE,STDOUT
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
URL = os.environ.get('DHUS_URL')
if URL is None:
    URL = 'https://scihub.copernicus.eu/dhus'
WAIT_TIME = 300
ONLINE_CHECK_TIME = 300
MAX_RETRY = 100
N_REQUEST = 3

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-u','--user',default=None,help='Username (or environment variable DHUS_USER is set)')
parser.add_option('-p','--password',default=None,help='Password (or environment variable DHUS_PASSWORD is set)')
parser.add_option('-U','--url',default=URL,help='API URL (or environment variable DHUS_URL is set)')
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
parser.add_option('-W','--wait_time',default=WAIT_TIME,type='int',help='Wait time to download data in sec (%default)')
parser.add_option('-O','--online_check_time',default=ONLINE_CHECK_TIME,type='int',help='Wait time to check online data in sec (%default)')
parser.add_option('-M','--max_retry',default=MAX_RETRY,type='int',help='Maximum number of retries to download data (%default)')
parser.add_option('-R','--n_request',default=N_REQUEST,type='int',help='Number of requests in advance (%default)')
parser.add_option('-d','--download',default=False,action='store_true',help='Download all results of the query. (%default)')
parser.add_option('-Y','--sort_year',default=False,action='store_true',help='Sort files by year. (%default)')
parser.add_option('-C','--no_checksum',default=False,action='store_true',help='Do NOT verify the downloaded files\' integrity by checking its MD5 checksum. (%default)')
parser.add_option('-f','--footprints',default=False,action='store_true',help='Create a geojson file search_footprints.geojson with footprints and metadata of the returned products. (%default)')
parser.add_option('-v','--version',default=False,action='store_true',help='Show the version and exit. (%default)')
parser.add_option('-Q','--quiet',default=False,action='store_true',help='Quiet mode (%default)')
(opts,args) = parser.parse_args()

def query_data(uuid):
    command = 'wget'
    command += ' --no-check-certificate'
    if opts.user is not None:
        command += ' --user {}'.format(opts.user)
    if opts.password is not None:
        command += ' --password {}'.format(opts.password)
    command += ' --output-document -'
    command += ' "'+opts.url+'/odata/v1/Products(\'{}\')"'.format(uuid)
    out = check_output(command,shell=True,stderr=PIPE).decode()
    root = ET.fromstring(out)
    child = None
    for value in root:
        prefix,has_namespace,postfix = value.tag.rpartition('}')
        if postfix.lower() == 'properties':
            child = value
            break
    if child is None:
        raise ValueError('Error, faild in finding properies >>> {}'.format(uuid))
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
    md5 = tags['Checksum']['Value']
    return name,size,stat,md5

def download_data(uuid,dst):
    dnam = os.path.dirname(dst)
    if dnam: # != ''
        if not os.path.exists(dnam):
            os.makedirs(dnam)
        if not os.path.isdir(dnam):
            raise IOError('Error, no such directory >>> '+dnam)
    command = 'wget'
    #command += ' --content-disposition'
    command += ' --continue'
    if opts.user is not None:
        command += ' --user {}'.format(opts.user)
    if opts.password is not None:
        command += ' --password {}'.format(opts.password)
    if opts.quiet:
        command += ' --quiet'
    command += ' --output-document '+dst
    command += ' "'+opts.url+'/odata/v1/Products(\'{}\')/\$value"'.format(uuid)
    try:
        call(command,shell=True)
        ret = 0
    except Exception:
        ret = -1
    return ret

# Query products
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
#if opts.no_checksum:
#    command += ' --no-checksum'
if opts.footprints:
    command += ' --footprints'
if opts.version:
    command += ' --version'
sys.stderr.write(command+'\n')
sys.stderr.flush()
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

# Get details of products
names = []
sizes = []
stats = []
md5s = []
for i,uuid in enumerate(uuids):
    name,size,stat,md5 = query_data(uuid)
    names.append(name)
    sizes.append(size)
    stats.append(stat)
    md5s.append(md5)
    sys.stderr.write('{:4d} {:40s} {:70s} {:10d} {:7s}\n'.format(i+1,uuid,name,size,'Online' if stat else 'Offline'))
    sys.stderr.flush()

if opts.download:
    # Make download list
    fnams = []
    download_list = []
    path = '.' if opts.path is None else opts.path
    for i in range(len(uuids)):
        # Check data availability
        if opts.sort_year:
            m = re.search('^\S+[^\d]_('+'\d'*4+')'+'\d'*4+'T'+'\d'*6+'_',names[i])
            if not m:
                raise ValueError('Error in file name >>> '+names[i])
            year = m.group(1)
            dnam = os.path.join(path,year)
        else:
            dnam = path
        fnam = os.path.join(dnam,names[i]+'.zip')
        gnam = fnam+'.incomplete'
        fnams.append(fnam)
        flag = False
        if os.path.exists(fnam) and os.path.getsize(fnam) == sizes[i]:
            flag = True
        if os.path.exists(gnam) and os.path.getsize(gnam) == sizes[i]:
            flag = True
        if not flag:
            download_list.append((uuids[i],fnam))
    # Request data in advance
    for i in range(min(opts.n_request,len(download_list))):
        uuid = download_list[i][0]
        fnam = download_list[i][1]
        gnam = fnam+'.incomplete'
        name,size,stat,md5 = query_data(uuid)
        if not stat:
            download_data(uuid,gnam)
    for i in range(len(uuids)):
        # Request data in advance
        if fnams[i] in download_list:
            j = download_list.index(fnams[i])+opts.n_request
            if j < len(download_list):
                uuid = download_list[j][0]
                fnam = download_list[j][1]
                gnam = fnam+'.incomplete'
                name,size,stat,md5 = query_data(uuid)
                if not stat:
                    download_data(uuid,gnam)
        uuid = uuids[i]
        fnam = fnams[i]
        gnam = fnam+'.incomplete'
        # Rename if gnam with expected size exists
        if os.path.exists(gnam):
            gsiz = os.path.getsize(gnam)
            if gsiz == sizes[i]:
                os.rename(gnam,fnam)
        # Skip if fnam with expected size exists
        if os.path.exists(fnam):
            fsiz = os.path.getsize(fnam)
            if fsiz == sizes[i]:
                if not opts.no_checksum:
                    with open(fnam,'rb') as fp:
                        md5 = hashlib.md5(fp.read()).hexdigest()
                    if md5 != md5s[i]:
                        sys.stderr.write('Warning, md5={}, md5s[{}]={} >>> {}\n'.format(md5,i,md5s[i],fnam))
                        sys.stderr.flush()
                sys.stderr.write('###### Successfully downloaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
                if os.path.exists(gnam):
                    os.remove(gnam)
                continue
        # Wait online
        ntry = 0
        while True:
            name,size,stat,md5 = query_data(uuid)
            if stat: # Online
                break
            sys.stderr.write('Offline. Wait for {} sec >>> {}\n'.format(opts.online_check_time,fnam))
            sys.stderr.flush()
            # Request data
            if ntry%opts.max_retry == 0:
                download_data(uuid,gnam)
            time.sleep(opts.online_check_time)
            ntry += 1
            continue
        # Download data
        for ntry in range(opts.max_retry): # loop to download 1 file
            download_data(uuid,gnam)
            # Rename if gnam with expected size exists
            if os.path.exists(gnam):
                gsiz = os.path.getsize(gnam)
                if gsiz == sizes[i]:
                    os.rename(gnam,fnam)
                elif gsiz > sizes[i]:
                    os.remove(gnam)
            # Exit if fnam with expected size exists
            if os.path.exists(fnam):
                fsiz = os.path.getsize(fnam)
                if fsiz == sizes[i]:
                    if not opts.no_checksum:
                        with open(fnam,'rb') as fp:
                            md5 = hashlib.md5(fp.read()).hexdigest()
                        if md5 != md5s[i]:
                            sys.stderr.write('Warning, md5={}, md5s[{}]={} >>> {}\n'.format(md5,i,md5s[i],fnam))
                            sys.stderr.flush()
                    sys.stderr.write('###### Successfully downloaded >>> {}\n'.format(fnam))
                    sys.stderr.flush()
                    if opts.log is not None:
                        with open(opts.log,'a') as fp:
                            fp.write(fnam+'\n')
                    if os.path.exists(gnam):
                        os.remove(gnam)
                    break
            sys.stderr.write('Wait for {} sec\n'.format(opts.wait_time))
            sys.stderr.flush()
            time.sleep(opts.wait_time)