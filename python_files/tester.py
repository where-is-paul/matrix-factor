#! /usr/bin/python
import re
import urllib
import sys
import os
import time
import tarfile
import shutil
import subprocess

def download_file(download_url, file_name):
    """Given urls, downloads files to current directory.
    """
    
    print 'Retrieving... ', download_url
    try:
        print file_name
        urllib.urlretrieve(download_url, file_name)
        print 'Done.'
        time.sleep(1)
        return True
    except IOError:
        print 'File not found.'
        return False
        
            
def find_files(text, file_type):
    file_list = re.findall(r'A HREF="([^<>]*?\.{})"'.format(file_type), text, re.IGNORECASE)
    return file_list

def main():
    args = sys.argv[1:]

    if not args:
        print 'usage: file_name'
        sys.exit(1)
        
    file_name = args[0]
    print '\nScanning for links from: ' + file_name + '\n'
    
    ufile = open(file_name)
    text = ufile.read()
    
    dlist = find_files(text, 'tar.gz')
    preprocess = lambda x: "MM" in x
    dlist = [link for link in dlist if preprocess(link)]
    print 'Files being downloaded:', '\n' + '\n'.join(dlist) + '\n'

    for download_url in dlist:
        file_name = os.path.split(download_url)[1]
        if download_file(download_url, file_name):
            tar = tarfile.open(file_name)
            tar.extractall()
            tar.close()

            file_no_ext = file_name.split('.')[0]
            mtx_path = '{}/{}.mtx'.format(file_no_ext, file_no_ext)
            output = subprocess.check_output(["../ldl_driver {} 1.0 0.001 -n -n".format(mtx_path)], shell=True)
            print output

            print 'Removing files... '
            shutil.rmtree('{}'.format(file_no_ext))
            os.remove(file_name)
    
if __name__ == '__main__':
    main()
