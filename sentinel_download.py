#!/usr/bin/env python
import os
import sys
from subprocess import call
from optparse import OptionParser,IndentedHelpFormatter

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
parser.add_option('-d','--download',default=False,action='store_true',help='Download all results of the query. (%default)')
parser.add_option('-n','--no_checksum',default=False,action='store_true',help='Do NOT verify the downloaded files\' integrity by checking its MD5 checksum. (%default)')
parser.add_option('-f','--footprints',default=False,action='store_true',help='Create a geojson file search_footprints.geojson with footprints and metadata of the returned products. (%default)')
parser.add_option('-v','--version',default=False,action='store_true',help='Show the version and exit. (%default)')
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
if opts.download:
    command += ' --download'
if opts.no_checksum:
    command += ' --no_checksum'
if opts.footprints:
    command += ' --footprints'
if opts.version:
    command += ' --version'
print(command)
