import os
from os.path import isdir , isfile , join

import numpy as np
import astropy as ap

from astropy.io import fits



if __name__ == '__main__':
	
	# getting all directories in current working dir
	wd = './'
	dir_names = [wd + f + '/' for f in os.listdir (wd) if isdir(join(wd,f))]

	star_params = np.loadtxt ('star_parameters.csv', skiprows=1 , dtype=str, delimiter=',')
	# print(star_params[:,0])
	cmd_list = []
	for d in dir_names : 
		if '00_m_16' in d or 'small_asteroid_lightcurve_CFHT_data' in d: continue

		file_names = [d + f for f in os.listdir(d) if isfile(join(d,f))]
		

		for f in file_names:
			# only working on new files around asteroid chip
			if not 'on' in f: continue
			if not ('EL157' in f ): continue

			try:
				file = fits.open (f)
				# print( f'{f} is a fits file')
			except Exception as e:
				# print( f'{f} is not a fits file' )
				continue

			L = 0 
			a = 0
			oid = f.split('_')[1]
			# print(oid)
			for i in range(len(star_params)):
				if oid in star_params[i,0]: 
					L = star_params[i,1]
					a = star_params[i,2]
					break

			command = f'\n{f[:]}'
			print(command)
			cmd_list.append(command)
	
	# for i in cmd_list: print(i)
	# output = open('ast_files.txt', 'w+')
	output = open('fits_files_all_chips.txt', 'w+')
	output.writelines( cmd_list )


