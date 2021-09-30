#!/usr/bin/env python
import os
import sys
import shutil
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
parser.add_option('-S','--srcdir',default=None,help='Source directory (%default)')
parser.add_option('-D','--dstdir',default=None,help='Destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()

opts.srcdir = os.path.abspath(opts.srcdir)
topdir = os.getcwd()
os.chdir(opts.drvdir)

folders = {}

def make_folder(path):
    global folders
    parent = os.path.dirname(path)
    target = os.path.basename(path)
    if parent == '':
        l = drive.ListFile({'q': '"root" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(target)}).GetList()
        n_list = len(l)
        if n_list != 1:
            os.chdir(topdir)
            raise ValueError('Error in finding folder, n_list={} >>> {}'.format(n_list,path))
        folders.update({path:l[0]})
        return 0
    elif not parent in folders:
        os.chdir(topdir)
        raise IOError('Error, no such folder >>> '+parent)
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType = "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list > 1:
        os.chdir(topdir)
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

def copy_file(fnam,gnam):
    parent = os.path.dirname(gnam)
    target = os.path.basename(gnam)
    if parent == '':
        os.chdir(topdir)
        raise ValueError('Error in parent >>> {}'.format(gnam))
    elif not parent in folders:
        os.chdir(topdir)
        raise IOError('Error, no such folder >>> '+parent)
    l = drive.ListFile({'q': '"{}" in parents and trashed = false and mimeType != "application/vnd.google-apps.folder" and title = "{}"'.format(folders[parent]['id'],target)}).GetList()
    n_list = len(l)
    if n_list > 1:
        os.chdir(topdir)
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
            return 0
    f = drive.CreateFile({'parents':[{'id':folders[parent]['id']}]})
    f.SetContentFile(fnam)
    f['title'] = target
    f.Upload()
    return 0

gauth = GoogleAuth()
gauth.LocalWebserverAuth()
drive = GoogleDrive(gauth)

for root,ds,fs in os.walk(opts.srcdir):
    curdir = os.path.relpath(root,opts.srcdir)
    if opts.verbose:
        sys.stderr.write('#####################\n')
        sys.stderr.write(curdir+'\n')
        sys.stderr.write(str(ds)+'\n')
        sys.stderr.write(str(fs)+'\n')
        sys.stderr.flush()
    if curdir == os.curdir:
        srcdir = opts.srcdir
        dstdir = opts.dstdir
    else:
        srcdir = os.path.join(opts.srcdir,curdir)
        dstdir = os.path.join(opts.dstdir,curdir)
    #print(srcdir,'-----',dstdir)
    if not dstdir in folders:
        if make_folder(dstdir) != 0:
            raise IOError('Error, faild in making folder >>> '+dstdir)
    for f in fs:
        fnam = os.path.join(srcdir,f)
        gnam = os.path.join(dstdir,f)
        copy_file(fnam,gnam)
os.chdir(topdir)
