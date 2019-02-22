import progressbar
import os

from nco import Nco


nco = Nco()

def get_files(mypath):
    """Returns a list of files that have filenames starting with "firstnine" as list of lists [fullpath, filename]."""
    fileinfo = []
    for root, dirs, files in os.walk(mypath, topdown = False):
        for name in files:
            if name[-2:] == 'nc':
                fileinfo.append([os.path.join(root,name),name])
    return (fileinfo)

def transpose_netcdf(path = '/data/marble/sandbox/jsilberstein/', target = '/data/marble/sandbox/jsilberstein/transposed/'):
	for file in progressbar.progressbar(get_files(path)):
		nco.ncpdq(input = file[0], output = '{}{}'.format(target,file[1]), arrange = ['latitude', 'longitude','time'])