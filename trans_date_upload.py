#!/usr/bin/env python
import os
import sys
import re
from datetime import datetime
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from optparse import OptionParser,IndentedHelpFormatter

# Defaults
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('USERPROFILE')
DRVDIR = os.path.join(HOME,'Work','SATREPS','IPB_Satreps')
SITE = 'Cihea'
LEVEL = 'test'
VERSION = 'v1.0'

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('--site',default=SITE,help='Site name (%default)')
parser.add_option('--level',default=LEVEL,help='Analysis level, test/final/preliminary (%default)')
parser.add_option('--version',default=VERSION,help='Product version (%default)')
parser.add_option('--date',default=None,help='Date in the format YYYYMMDD (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
if opts.date is None:
    opts.date = datetime.now().strftime('%Y%m%d')
fnams = args

topdir = os.getcwd()
os.chdir(opts.drvdir)

gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

# Get Spatial-Information folder
l = drive.ListFile({'q': '"root" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "Spatial-Information"'}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding Spatial-Information folder')
folder_spatial_information = l[0]
# Get Transplanting_date folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "Transplanting_date"'.format(folder_spatial_information['id'])}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding Transplanting_date folder')
folder_transplanting_date = l[0]
# Get site folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_transplanting_date['id'],opts.site)}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding {} folder'.format(opts.site))
folder_site = l[0]
# Get level folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_site['id'],opts.level.lower())}).GetList()
if len(l) != 1:
    raise ValueError('Error in finding {} folder'.format(opts.level))
folder_level = l[0]
# Get version folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_level['id'],opts.version)}).GetList()
if len(l) != 1:
    folder_version = drive.CreateFile({'parents':[{'id':folder_level['id']}],'mimeType':'application/vnd.google-apps.folder','title':opts.version})
    folder_version.Upload()
else:
    folder_version = l[0]
# Get year folder
dstr_year = datetime.strptime(opts.date,'%Y%m%d').strftime('%Y')
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_version['id'],dstr_year)}).GetList()
if len(l) != 1:
    folder_year = drive.CreateFile({'parents':[{'id':folder_version['id']}],'mimeType':'application/vnd.google-apps.folder','title':dstr_year})
    folder_year.Upload()
else:
    folder_year = l[0]
# Get date folder
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title contains "{}"'.format(folder_year['id'],opts.date)}).GetList()
if len(l) != 1:
    folder_date = drive.CreateFile({'parents':[{'id':folder_year['id']}],'mimeType':'application/vnd.google-apps.folder','title':opts.date})
    folder_date.Upload()
else:
    folder_date = l[0]
l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType != "application/vnd.google-apps.folder"'.format(folder_date['id'])}).GetList()

for input_fnam in fnams:
    upload_fnam = os.path.basename(input_fnam)
    flag = True
    for f in l:
        if f['title'].upper() == upload_fnam.upper():
            if opts.overwrite:
                sys.stderr.write('File exists, delete   >>> '+f['title']+'\n')
                f.Delete()
            else:
                sys.stderr.write('File exists, skip     >>> '+f['title']+'\n')
                flag = False # no need to upload
                break
    if flag: # upload file
        f = drive.CreateFile({'parents':[{'id':folder_date['id']}]})
        f.SetContentFile(input_fnam)
        f['title'] = upload_fnam
        f.Upload()
        sys.stderr.write('Successfully uploaded >>> '+upload_fnam+'\n')

os.chdir(topdir)
