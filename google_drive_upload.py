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
parser.add_option('-S','--srcdir',default=None,help='Source directory (%default)')
parser.add_option('-s','--subdir',default=None,action='append',help='Sub directory (%default)')
parser.add_option('-D','--dstdir',default=None,help='Destination directory (%default)')
parser.add_option('-L','--locdir',default=None,help='Local destination directory (%default)')
parser.add_option('--drvdir',default=DRVDIR,help='GoogleDrive directory (%default)')
parser.add_option('-K','--keep_folder',default=None,action='append',help='Directory to keep (%default)')
parser.add_option('-I','--ignore_file',default=None,action='append',help='File to ignore (%default)')
parser.add_option('-v','--verbose',default=False,action='store_true',help='Verbose mode (%default)')
parser.add_option('-d','--debug',default=False,action='store_true',help='Debug mode (%default)')
parser.add_option('--overwrite',default=False,action='store_true',help='Overwrite mode (%default)')
(opts,args) = parser.parse_args()
if opts.srcdir is None or opts.subdir is None or opts.dstdir is None or opts.locdir is None:
    raise ValueError('Error, srcdir={}, subdir={}, dstdir={}, locdir={}'.format(opts.srcdir,opts.subdir,opts.dstdir,opts.locdir))
keep_folder = []
keep_folder_lower = []
if opts.keep_folder is not None:
    for f in opts.keep_folder:
        for p in glob(os.path.normpath(os.path.join(opts.srcdir,f))):
            if os.path.isdir(p):
                keep_folder.append(p)
                keep_folder_lower.append(p.lower())
ignore_file = []
ignore_file_lower = []
if opts.ignore_file is not None:
    for f in opts.ignore_file:
        for p in glob(os.path.normpath(os.path.join(opts.srcdir,f))):
            if not os.path.isdir(p):
                ignore_file.append(p)
                ignore_file_lower.append(p.lower())
if opts.verbose:
    if len(keep_folder) > 0:
        sys.stderr.write('keep_folder:\n')
        for f in keep_folder:
            sys.stderr.write(f+'\n')
        sys.stderr.flush()
    if len(ignore_file) > 0:
        sys.stderr.write('ignore_file:\n')
        for f in ignore_file:
            sys.stderr.write(f+'\n')
        sys.stderr.flush()

opts.srcdir = os.path.abspath(opts.srcdir)
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

for subdir in opts.subdir:
    make_folders(os.path.join(opts.dstdir,subdir))
    for root,ds,fs in os.walk(os.path.join(opts.srcdir,subdir)):
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
            locdir = opts.locdir
        else:
            srcdir = os.path.join(opts.srcdir,curdir)
            dstdir = os.path.join(opts.dstdir,curdir)
            locdir = os.path.join(opts.locdir,curdir)
        #print(srcdir,'-----',dstdir)
        if not dstdir in folders:
            if make_folder(dstdir) != 0:
                raise IOError('Error, faild in making folder >>> '+dstdir)
        if not os.path.exists(locdir):
            os.makedirs(locdir)
        if not os.path.isdir(locdir):
            raise IOError('Error, no such folder >>> '+locdir)
        for f in fs:
            fnam = os.path.join(srcdir,f)
            gnam = os.path.join(dstdir,f)
            lnam = os.path.join(locdir,f)
            if fnam.lower() in ignore_file_lower:
                continue
            if upload_and_check_file(fnam,gnam) == 0:
                shutil.move(fnam,lnam)
                if opts.debug and not os.path.exists(fnam) and os.path.exists(lnam):
                    sys.stderr.write('Moved {} to {}\n'.format(fnam,lnam))
                    sys.stderr.flush()
        if len(os.listdir(srcdir)) == 0:
            if not srcdir.lower() in keep_folder_lower:
                os.rmdir(srcdir)
                if opts.debug and not os.path.exists(srcdir):
                    sys.stderr.write('Removed {}\n'.format(srcdir))
                    sys.stderr.flush()
