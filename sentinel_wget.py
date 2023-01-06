#!/usr/bin/env python
import os
import psutil
import sys
import re
import time
import hashlib
import numpy as np
import xml.etree.ElementTree as ET
from collections import OrderedDict
from datetime import datetime
from subprocess import call,check_output,PIPE,STDOUT
from optparse import OptionParser,IndentedHelpFormatter

# Constants
KB = 1024
MB = KB*1024
GB = MB*1024

# Defaults
URL = os.environ.get('DHUS_URL')
if URL is None:
    URL = 'https://scihub.copernicus.eu/dhus'
DOTBYTES = '10M'
WAIT_TIME = 300
QUERY_WAIT_TIME = 60
ONLINE_CHECK_TIME = 600
CLEANUP_TIME = 172800 # 2 days
MAX_RETRY = 100
N_REQUEST = 20

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.add_option('-u','--user',default=None,help='Username (or environment variable DHUS_USER is set)')
parser.add_option('-p','--password',default=None,help='Password (or environment variable DHUS_PASSWORD is set)')
parser.add_option('-U','--url',default=URL,help='API URL (or environment variable DHUS_URL is set)')
parser.add_option('-n','--netrc',default=None,help='Netrc file name (%default)')
parser.add_option('-m','--machine',default=None,help='Machine name (%default)')
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
parser.add_option('-B','--dotbytes',default=DOTBYTES,help='Dot size of progress indicator (%default)')
parser.add_option('-L','--log',default=None,help='Set the log file name.')
parser.add_option('-P','--path',default=None,help='Set the path where the files will be saved.')
parser.add_option('-q','--query',default=None,help='Extra search keywords you want to use in the query. Separate keywords with comma. Example: \'producttype=GRD,polarisationmode=HH\'.')
parser.add_option('-W','--wait_time',default=WAIT_TIME,type='int',help='Wait time to download data in sec (%default)')
parser.add_option('-w','--query_wait_time',default=QUERY_WAIT_TIME,type='int',help='Wait time to query data in sec (%default)')
parser.add_option('-O','--online_check_time',default=ONLINE_CHECK_TIME,type='int',help='Wait time to check online data in sec (%default)')
parser.add_option('-T','--cleanup_time',default=CLEANUP_TIME,type='int',help='Time to cleanup incomplete files in sec (%default)')
parser.add_option('-M','--max_retry',default=MAX_RETRY,type='int',help='Maximum number of retries to download data (%default)')
parser.add_option('-R','--n_request',default=N_REQUEST,type='int',help='Number of requests in advance (%default)')
parser.add_option('-D','--dosmode',default=None,action='store_true',help='DOS mode (%default)')
parser.add_option('-d','--download',default=False,action='store_true',help='Download all results of the query. (%default)')
parser.add_option('-Y','--sort_year',default=False,action='store_true',help='Sort files by year. (%default)')
parser.add_option('-C','--no_checksum',default=False,action='store_true',help='Do NOT verify the downloaded files\' integrity by checking its MD5 checksum. (%default)')
parser.add_option('-f','--footprints',default=False,action='store_true',help='Create a geojson file search_footprints.geojson with footprints and metadata of the returned products. (%default)')
parser.add_option('-v','--version',default=False,action='store_true',help='Show the version and exit. (%default)')
parser.add_option('-Q','--quiet',default=False,action='store_true',help='Quiet mode (%default)')
(opts,args) = parser.parse_args()
if opts.netrc is not None:
    server = None
    username = None
    password = None
    flag = False
    with open(opts.netrc,'r') as fp:
        for line in fp:
            m = re.search('machine\s+(\S+)',line)
            if m:
                if re.search(opts.machine,m.group(1)):
                    flag = True
                    server = m.group(1)
                else:
                    flag = False
                continue
            m = re.search('login\s+(\S+)',line)
            if m:
                if flag:
                    username = m.group(1)
                continue
            m = re.search('password\s+(\S+)',line)
            if m:
                if flag:
                    password = m.group(1)
                continue
    if server is None or username is None or password is None:
        raise ValueError('Error, server={}, username={}, password={}'.format(server,username,password))
    opts.user = username
    opts.password = password

def clean_up():
    if not np.all(np.array([len(fnams),len(sizes),len(md5s),len(size_values),len(size_errors)]) == len(uuids)):
        raise ValueError('Error, different size.')
    tcur = time.time()
    for i in range(len(uuids)):
        fnam = fnams[i]
        gnam = fnam+'.incomplete'
        if os.path.exists(gnam):
            gsiz = os.path.getsize(gnam)
            gtim = os.path.getmtime(gnam)
            if (gsiz == sizes[i]) or (np.abs(gsiz-size_values[i]) <= size_errors[i]):
                pass # keep
            elif tcur-gtim > opts.cleanup_time:
                os.remove(gnam)
    return

def make_list():
    if not np.all(np.array([len(fnams),len(sizes),len(md5s),len(size_values),len(size_errors)]) == len(uuids)):
        raise ValueError('Error, different size.')
    d_list = []
    for i in range(len(uuids)):
        uuid = uuids[i]
        fnam = fnams[i]
        gnam = fnam+'.incomplete'
        flag = False
        if os.path.exists(fnam):
            fsiz = os.path.getsize(fnam)
            if (fsiz == sizes[i]) or (np.abs(fsiz-size_values[i]) <= size_errors[i]):
                flag = True
        if os.path.exists(gnam):
            gsiz = os.path.getsize(gnam)
            if (gsiz == sizes[i]) or (np.abs(gsiz-size_values[i]) <= size_errors[i]):
                flag = True
        if not flag:
            d_list.append((uuid,fnam))
    return d_list

def check_data(i,e_list):
    if not np.all(np.array([len(fnams),len(sizes),len(md5s),len(size_values),len(size_errors)]) == len(uuids)):
        raise ValueError('Error, different size.')
    fnam = fnams[i]
    gnam = fnam+'.incomplete'
    # Rename if gnam with expected size exists
    if os.path.exists(gnam):
        gsiz = os.path.getsize(gnam)
        if (gsiz == sizes[i]) or (np.abs(gsiz-size_values[i]) <= size_errors[i]):
            os.rename(gnam,fnam)
        elif gsiz > max(sizes[i],size_values[i]+size_errors[i]):
            os.remove(gnam)
    if os.path.exists(fnam):
        fsiz = os.path.getsize(fnam)
        if (fsiz == sizes[i]) or (np.abs(fsiz-size_values[i]) <= size_errors[i]):
            if not opts.no_checksum:
                with open(fnam,'rb') as fp:
                    md5 = hashlib.md5(fp.read()).hexdigest()
                if md5.upper() != md5s[i].upper():
                    sys.stderr.write('Warning, md5={}, md5s[{}]={} >>> {}\n'.format(md5,i,md5s[i],fnam))
                    sys.stderr.flush()
            sys.stderr.write('###### File exists >>> {}\n'.format(fnam))
            sys.stderr.flush()
            if opts.log is not None and not fnam in e_list:
                with open(opts.log,'a') as fp:
                    fp.write(fnam+'\n')
            if os.path.exists(gnam):
                os.remove(gnam)
            return 0
    return -1

def query_data(uuid):
    command = 'wget'
    command += ' --no-check-certificate'
    if opts.user is not None:
        command += ' --user "{}"'.format(opts.user)
    if opts.password is not None:
        command += ' --password "{}"'.format(opts.password)
    command += ' --output-document -'
    command += ' "'+opts.url+'/odata/v1/Products(\'{}\')"'.format(uuid)
    while True:
        try:
            out = check_output(command,shell=True,stderr=PIPE).decode()
            break
        except Exception as inst:
            sys.stderr.write(str(inst)+'\n')
            sys.stderr.write('Query failed. Wait for {} sec.\n'.format(opts.query_wait_time))
            sys.stderr.flush()
            time.sleep(opts.query_wait_time)
    root = ET.fromstring(out)
    child = None
    for value in root:
        prefix,has_namespace,postfix = value.tag.rpartition('}')
        if postfix.lower() == 'properties':
            child = value
            break
    if child is None:
        raise ValueError('Error, failed in finding properies >>> {}'.format(uuid))
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

def download_data(uuid,fnam):
    if not np.all(np.array([len(fnams),len(sizes),len(md5s),len(size_values),len(size_errors)]) == len(uuids)):
        raise ValueError('Error, different size.')
    elif not fnam in fnams:
        raise ValueError('Error, not in fnams >>> '+fnam)
    else:
        indx = fnams.index(fnam)
    gnam = fnam+'.incomplete'
    dnam = os.path.dirname(fnam)
    if dnam: # != ''
        if not os.path.exists(dnam):
            os.makedirs(dnam)
        if not os.path.isdir(dnam):
            raise IOError('Error, no such directory >>> '+dnam)
    command = 'wget'
    #command += ' --content-disposition'
    command += ' --no-check-certificate'
    command += ' --continue'
    command += ' --progress=dot'
    command += ' --execute dotbytes={}'.format(opts.dotbytes)
    if opts.user is not None:
        command += ' --user "{}"'.format(opts.user)
    if opts.password is not None:
        command += ' --password "{}"'.format(opts.password)
    if opts.quiet:
        command += ' --quiet'
    command += ' --output-document '+gnam
    if opts.dosmode:
        command += ' "'+opts.url+'/odata/v1/Products(\'{}\')/$value"'.format(uuid)
    else:
        command += ' "'+opts.url+'/odata/v1/Products(\'{}\')/\$value"'.format(uuid)
    try:
        sys.stderr.write('###### Download request for {}\n'.format(fnam))
        sys.stderr.flush()
        ret = call(command,shell=True)
        if ret != 0:
            raise ValueError('Error, command returned non-zero exit status >>> {}'.format(ret))
        ret = 1
        # Rename if gnam with expected size exists
        if os.path.exists(gnam):
            gsiz = os.path.getsize(gnam)
            if (gsiz == sizes[indx]) or (np.abs(gsiz-size_values[indx]) <= size_errors[indx]):
                os.rename(gnam,fnam)
            elif gsiz > max(sizes[indx],size_values[indx]+size_errors[indx]):
                os.remove(gnam)
        if os.path.exists(fnam):
            fsiz = os.path.getsize(fnam)
            if (fsiz == sizes[indx]) or (np.abs(fsiz-size_values[indx]) <= size_errors[indx]):
                if not opts.no_checksum:
                    with open(fnam,'rb') as fp:
                        md5 = hashlib.md5(fp.read()).hexdigest()
                    if md5.upper() != md5s[indx].upper():
                        sys.stderr.write('Warning, md5={}, md5s[{}]={} >>> {}\n'.format(md5,indx,md5s[indx],fnam))
                        sys.stderr.flush()
                sys.stderr.write('###### Successfully downloaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
                ret = 0
    except Exception as inst:
        sys.stderr.write(str(inst)+'\n')
        sys.stderr.flush()
        ret = -1
    return ret

def request_data(d_list):
    for i in range(min(opts.n_request,len(d_list))):
        uuid = d_list[i][0]
        fnam = d_list[i][1]
        gnam = fnam+'.incomplete'
        if os.path.exists(gnam):
            continue
        name,size,stat,md5 = query_data(uuid)
        if not stat:
            ret = download_data(uuid,fnam)
            if ret < 0:
                if os.path.exists(gnam):
                    gsiz = os.path.getsize(gnam)
                    if gsiz == 0:
                        os.remove(gnam)
    return

def download_next_data(d_list):
    ret = -1
    for i in range(len(d_list)):
        uuid = d_list[i][0]
        fnam = d_list[i][1]
        name,size,stat,md5 = query_data(uuid)
        if stat:
            ret = download_data(uuid,fnam)
            if ret == 0:
                break
    return ret

if opts.dosmode is None:
    out = check_output('echo \$',shell=True).decode().strip()
    if len(out) == 2: # '\$'
        opts.dosmode = True
    else: # '$'
        opts.dosmode = False

# Query products
command = 'sentinelsat'
if opts.user is not None:
    command += ' --user "{}"'.format(opts.user)
if opts.password is not None:
    command += ' --password "{}"'.format(opts.password)
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
while True:
    try:
        out = check_output(command,shell=True,stderr=STDOUT).decode()
        break
    except Exception as inst:
        sys.stderr.write(str(inst)+'\n')
        sys.stderr.write('Search failed. Wait for {} sec.\n'.format(opts.query_wait_time))
        sys.stderr.flush()
        time.sleep(opts.query_wait_time)
sys.stderr.write(out+'\n')
sys.stderr.flush()
if opts.version:
    sys.exit()
uuids = []
size_labels = []
size_values = []
size_errors = []
for line in out.splitlines():
    # Product bb7a7783-f91b-4a35-907d-a6ddb807da73 - Date: 2019-07-09T02:55:59.024Z, Instrument: MSI, Mode: , Satellite: Sentinel-2, Size: 1.09 GB
    m = re.search('Product\s+(\S+)',line)
    if not m:
        continue
    uuid = m.group(1)
    m = re.search('Size\s*:\s+(\S+)\s+(\S+)',line)
    if not m:
        raise ValueError('Error in finding file size >>> '+line)
    label = m.group(1)+' '+m.group(2)
    value = float(m.group(1))
    unit = m.group(2).upper()
    if unit == 'KB':
        value *= KB
        error = KB//10
    elif unit == 'MB':
        value *= MB
        error = MB//10
    elif unit == 'GB':
        value *= GB
        error = GB//10
    else:
        raise ValueError('Error in size unit >>> '+line)
    uuids.append(uuid)
    size_labels.append(label)
    size_values.append(value)
    size_errors.append(error)

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
    if np.abs(size-size_values[i]) > size_errors[i]:
        size_flag = '!'
    else:
        size_flag = ' '
        size_values[i] = size
        size_errors[i] = 0
    sys.stderr.write('{:4d} {:40s} {:70s} {:10d} {} ({:>10s}) {:7s}\n'.format(i+1,uuid,name,size,size_flag,size_labels[i],'Online' if stat else 'Offline'))
    sys.stderr.flush()

if opts.download:
    # Make file list
    fnams = []
    path = '.' if opts.path is None else opts.path
    for i in range(len(uuids)):
        if opts.sort_year:
            m = re.search('^\S+[^\d]_('+'\d'*4+')'+'\d'*4+'T'+'\d'*6+'_',names[i])
            if not m:
                raise ValueError('Error in file name >>> '+names[i])
            year = m.group(1)
            dnam = os.path.join(path,year)
        else:
            dnam = path
        fnam = os.path.join(dnam,names[i]+'.zip')
        fnams.append(fnam)
    clean_up() # Remove old .incomplete files
    download_list = make_list() # Make download list
    existing_list = [fnam for fnam in fnams if not fnam in [l[1] for l in download_list]]
    request_data(download_list) # Request data in advance
    for i in range(len(uuids)):
        uuid = uuids[i]
        fnam = fnams[i]
        gnam = fnam+'.incomplete'
        # Skip if fnam or gnam with expected size exists
        if check_data(i,existing_list) == 0:
            continue
        # Wait online
        while True:
            name,size,stat,md5 = query_data(uuid)
            if stat: # Online
                break
            tpre = time.time()
            request_data([(uuid,fnam)]) # Request data
            if (download_next_data(download_list) == 0): # Download next data
                clean_up() # Remove old .incomplete files
                download_list = make_list() # Make download list
                request_data(download_list) # Request data in advance
            tdif = time.time()-tpre
            if tdif < opts.online_check_time:
                sys.stderr.write('Offline. Wait for {} sec >>> {}\n'.format(opts.online_check_time-tdif,fnam))
                sys.stderr.flush()
                time.sleep(opts.online_check_time-tdif)
            continue
        # Download data
        for ntry in range(opts.max_retry): # loop to download 1 file
            tpre = time.time()
            if (download_data(uuid,fnam) == 0):
                clean_up() # Remove old .incomplete files
                download_list = make_list() # Make download list
                request_data(download_list) # Request data in advance
            # Exit if fnam or gnam with expected size exists
            if check_data(i,existing_list) == 0:
                break
            tdif = time.time()-tpre
            if tdif < opts.wait_time:
                if (download_next_data(download_list) == 0): # Download next data
                    clean_up() # Remove old .incomplete files
                    download_list = make_list() # Make download list
                    request_data(download_list) # Request data in advance
                tdif = time.time()-tpre
                if tdif < opts.wait_time:
                    sys.stderr.write('Wait for {} sec\n'.format(opts.wait_time-tdif))
                    sys.stderr.flush()
                    time.sleep(opts.wait_time-tdif)
