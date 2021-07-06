#!/usr/bin/env python
import os
import sys
import re
from datetime import datetime
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
SITE = 'Test'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-S','--site',default=SITE,help='Target sites (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args

gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

# Get Spatial-Information folder
l = drive.ListFile({'q': '"root" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "Spatial-Information"'}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding Spatial-Information folder')
folder_spatial_information = l[0]
# Get SENTINEL-1 folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "SENTINEL-2"'.format(folder_spatial_information['id'])}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding SENTINEL-2 folder')
folder_sentinel_1 = l[0]
# Get L2A folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "L2A"'.format(folder_sentinel_1['id'])}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding L2A folder')
folder_grd = l[0]
# Get SITE folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_grd['id'],opts.site)}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding {} folder'.format(opts.site))
folder_site = l[0]

for input_fnam in fnams:
    # S2A_MSIL2A_20210104T030121_N0214_R032_T48MYT_20210104T062157.zip
    fnam = os.path.basename(input_fnam)
    unam = fnam.upper()
    bnam,enam = os.path.splitext(unam)
    if enam != '.ZIP':
        continue
    m = re.search('([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)',unam)
    if not m:
        continue
    if m.group(1) != 'S2A' and m.group(1) != 'S2B':
        continue
    if m.group(2) != 'MSIL2A':
        sys.stderr.write('Warning, skipping file >>> '+fnam+'\n')
        continue
    d1 = datetime.strptime(m.group(3),'%Y%m%dT%H%M%S')
    d2 = datetime.strptime(m.group(7),'%Y%m%dT%H%M%S')
    if d1.date() != d2.date():
        raise ValueError('Error, d1={}, d2={}'.format(m.group(3),m.group(7)))
    search_key = m.group(1)+'_'+m.group(2)+'_'+d1.strftime('%Y%m%d')
    upload_fnam = bnam+enam.lower()
    dstr_year = d1.strftime('%Y')
    # Get Year folder
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_site['id'],dstr_year)}).GetList()
    if len(l) != 1:
        folder_year = drive.CreateFile({'parents':[{'id':folder_site['id']}],'mimeType':'application/vnd.google-apps.folder','title':dstr_year})
        folder_year.Upload()
    else:
        folder_year = l[0]
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType != "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_year['id'],search_key)}).GetList()
    flag = True
    for f in l:
        if f['title'].upper() == unam:
            if opts.overwrite:
                sys.stderr.write('File exists, delete   >>> '+f['title']+'\n')
                f.Delete()
            else:
                sys.stderr.write('File exists, skip     >>> '+f['title']+'\n')
                flag = False # no need to upload
                break
        else:
            sys.stderr.write('Warning, different file for the same date >>> '+f['title']+'\n')
    if flag: # upload file
        f = drive.CreateFile({'parents':[{'id':folder_year['id']}]})
        f.SetContentFile(input_fnam)
        f['title'] = upload_fnam
        f.Upload()
        sys.stderr.write('Successfully uploaded >>> '+fnam+'\n')
