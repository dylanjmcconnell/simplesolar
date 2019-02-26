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
    fileinfo.sort()
    return (fileinfo)


def transpose_netcdf(path = '/data/marble/sandbox/jsilberstein/years/', target = '/data/marble/sandbox/jsilberstein/yearstransposed/'):
	"""This transposes the netcdf to make time the fastest reading demension. Note that if time is fastest reading demension concatenation takes OOM longer to run."""
	for file in progressbar.progressbar(get_files(path)):
		nco.ncpdq(input = file[0], output = '{}llt_{}'.format(target,file[1]), arrange = ['latitude', 'longitude','time'])


def convert_time_to_double(path = '/data/marble/sandbox/jsilberstein/', target = '/data/marble/sandbox/jsilberstein/timedouble/'):
	"""This converts the time unit to a double to ensure that the the concatenation offsets the time appropriately (time as integer causes date to restart (01/01/01 00:00:00) each month)."""
	for file in progressbar.progressbar(get_files(path)):
		nco.ncap2(options= ['-O','-s','time=time.double()'], input =  file[0], output = '{}dbl_tm_{}'.format(target,file[1])) 

def concat_files(path = '/data/marble/sandbox/jsilberstein/timedouble/', target = '/data/marble/sandbox/jsilberstein/timedoubleyears/'):
	"""This file concatenates multiple files along the time dimension. Recommend having time as the first dim in order for speed."""
	year_lists = []
	year = []
	files = get_files(path)
	for x in files:
		year.append(x[0])
		if x[0][-5:-3] == '12':
			year_lists.append(year)
			year = []
	for x in progressbar.progressbar(year_lists):
		nco.ncrcat(input = x, output = '{}radiation_{}.nc'.format(target, x[0][-9:-5]))

	return year_lists 
