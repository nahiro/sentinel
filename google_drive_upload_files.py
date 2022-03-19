#!/usr/bin/env python
import os
import sys
import shutil
import hashlib
import atexit
from glob import glob
from pydrive2.auth import GoogleAuth
from pydrive2.drive import GoogleDrive
from optparse import OptionParser,IndentedHelpFormatter

# Default values
HOME = os.environ.get('HOME')
if HOME is None:
    HOME = os.environ.get('HOMEPATH')
DRVDIR = os.path.join(HOME,'Work','SATREPS','IPB_Satreps')

# Read options
parser = OptionParser(formatter=IndentedHelpFormatter(max_help_position=200,width=200))
parser.set_usage('Usage: %prog list_of_input_file [options]')
parser.add_option('-D','--dstdir',default=None,help='Destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if len(args) < 1:
    parser.print_help()
    sys.exit(0)
fnams = args
if opts.dstdir is None:
    raise ValueError('Error, dstdir={}'.format(opts.dstdir))

topdir = os.getcwd()
os.chdir(opts.drvdir)
def clean_up():
    os.chdir(topdir)
atexit.register(clean_up)

folders = {}

def make_folder(path):
    global folders
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    if parent == '':
        l = drive.ListFile({'q': '"root" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(target)}).GetList()
        n_list = len(l)
        if n_list != 1:
            raise ValueError('Error in finding folder, n_list={} >>> {}'.format(n_list,path))
        folders.update({path:l[0]})
        return 0
    elif not parent in folders:
        raise IOError('Error, no such folder >>> '+parent)
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list > 1:
        raise ValueError('Error in finding folder, n_list={} >>> {}'.format(n_list,path))
    elif n_list == 1:
        if not path in folders:
            folders.update({path:l[0]})
        return 0
    else:
        folder = drive.CreateFile({'parents':[{'id':folders[parent]['id']}],'mimeType':'application/vnd.google-apps.folder','title':target})
        folder.Upload()
        folders.update({path:folder})
        return 0
    return -1

def make_folders(path):
    normalized_path = os.path.normpath(path)
    path_components = normalized_path.split(os.sep)
    dnam = ''
    for p in path_components:
        dnam = os.path.join(dnam,p)
        make_folder(dnam)

def upload_file(fnam,gnam):
    parent = os.path.dirname(gnam)
    target = os.path.basename(gnam)
    if parent == '':
        raise ValueError('Error in parent >>> {}'.format(gnam))
    elif not parent in folders:
        raise IOError('Error, no such folder >>> '+parent)
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType != "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list > 1:
        raise ValueError('Error, n_list={} >>> {}'.format(n_list,gnam))
    elif n_list == 1:
        f = l[0]
        if opts.overwrite:
            if opts.verbose:
                sys.stderr.write('File exists, delete >>> '+f['title']+'\n')
                sys.stderr.flush()
            f.Delete()
        else:
            if opts.verbose:
                sys.stderr.write('File exists, skip   >>> '+f['title']+'\n')
                sys.stderr.flush()
            return f['title'],int(f['fileSize']),f['modifiedDate'],f['md5Checksum']
    # Upload file
    f = drive.CreateFile({'parents':[{'id':folders[parent]['id']}]})
    f.SetContentFile(fnam)
    f['title'] = target
    f.Upload()
    # Check uploaded file
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType != "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list != 1:
        raise ValueError('Error, n_list={} >>> {}'.format(n_list,gnam))
    else:
        f = l[0]
        return f['title'],int(f['fileSize']),f['modifiedDate'],f['md5Checksum']

def upload_and_check_file(fnam,gnam):
    title = os.path.basename(gnam)
    size = os.path.getsize(fnam)
    title_dst,size_dst,mdate_dst,md5_dst = upload_file(fnam,gnam)
    if (title_dst.lower() == title.lower()) and (size_dst == size):
        with open(fnam,'rb') as fp:
            md5 = hashlib.md5(fp.read()).hexdigest()
        if (md5_dst.lower() == md5.lower()):
            if opts.verbose:
                sys.stderr.write('Successfully uploaded >>> {}\n'.format(fnam))
                sys.stderr.flush()
            return 0
        else:
            if opts.verbose:
                sys.stderr.write('Warning, title={} ({}), size={} ({}), md5={} ({})\n'.format(title_dst,title,size_dst,size,md5_dst,md5))
                sys.stderr.flush()
            return -1
    else:
        if opts.verbose:
            sys.stderr.write('Warning, title={} ({}), size={} ({})\n'.format(title_dst,title,size_dst,size))
            sys.stderr.flush()
        return -1

gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

# Upload file
make_folders(opts.dstdir)
for fnam in fnams:
    f = os.path.basename(fnam)
    gnam = os.path.join(opts.dstdir,f)
    upload_and_check_file(fnam,gnam)
