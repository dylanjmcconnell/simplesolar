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
    			fileinfo.append([os.path.join(root,name), name])
    return (fileinfo)


def transpose_netcdf(path = '/data/marble/sandbox/jsilberstein/', target = '/data/marble/sandbox/jsilberstein/transposed/'):
	for file in progressbar.progressbar(get_files(path)):
		nco.ncpdq(input = file[0], output = '{}{}'.format(target,file[1]), arrange = ['latitude', 'longitude','time'])

# def combine_netcdf(path = '/data/marble/sandbox/jsilberstein/transposed/'):
	# for year in files:
	# nco.ncrcat()

def convert_time_to_double(path = '/data/marble/sandbox/jsilberstein/', target = '/data/marble/sandbox/jsilberstein/timedouble/'):
	for file in progressbar.progressbar(get_files(path)):
		nco.ncap2(options= ['-O','-s','time=time.double()'], input =  file[0], output = '{}dbl_tm_{}'.format(target,file[1])) 

def concat_files(path = '/data/marble/sandbox/jsilberstein/timedouble/', target = '/data/marble/sandbox/jsilberstein/timedoubleyears/'):
	year_lists = []
	year = []
	files = get_files(path)
	print (files)
	for x in files:
		year.append(x[0])
		if x[1][-5:-3] == 12:
			year_lists.append(year)
			year = []
	# for x in progressbar.progressbar(year_lists):
	# 	nco.ncrcat(input = x, output = '{}radiation_{}.nc'.format(target, x[1][-9:-5]))

	return year_lists 
