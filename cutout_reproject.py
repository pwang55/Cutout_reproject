'''

Usage:

    In python script directory:
    $ python cutout_reproject.py path/masterlist.txt (4400)

    In stack files directory:
    $ python path_to_script/cutout_reproject.py masterlist.txt (4400)

This script takes masterlist.txt as inlut (which is a list of all the individual stacks from 1ne to 4sw).
It will find the master_stack_cluster.fits and make a cutout of it, then reproject all individual stacks to the master cutout.
The second argument 4400 is optional and it is the dimension of master cutout. Default is 4400 pixels.

'''
import numpy as np
from astropy.io import fits
from astropy import wcs
import sys
import time
import os.path
import cutout_reproject_fxs

np = 4400

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print(__doc__)
    exit()
else:
    listpathname = sys.argv[1]  # path/masterlist.txt
    if len(sys.argv) == 3:
        np = int(sys.argv[2])


with open(listpathname) as f:
    fname = f.readline().strip()        # get the first filename in the list
    clustername = fname.split('_')[1]   # get the cluster name from the first filename

listname = listpathname.split('/')[-1]  # masterlist.txt
listnamelen = len(listname)
listpath = listpathname[:-listnamelen]  # path of masterlist.txt that ends with "/"
filepath = listpath


# Make master cutout image
if os.path.isfile(listpath + 'master_stack_{}_cutout.fits'.format(clustername)):
    print('Master Stack Cutout already exist.')
else:
    cutout_reproject_fxs.master_cutout(path=listpath, filename='master_stack_{}.fits'.format(clustername), np=np)


# Reproject individual stack images to the master cutout image
with open(listpathname) as f:
    for l in f:
        filename = l.strip()
        if os.path.isfile('cutout_' + filepath + filename):
            print('Reproject cutout already exist:\t' + filename)
        else:
            cutout_reproject_fxs.reproject2mastercutout(path=filepath, filename=filename, master_cutout='master_stack_{}_cutout.fits'.format(clustername))


