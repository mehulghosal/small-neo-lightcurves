import numpy as np
import astropy as ap
from astropy.io import fits

import os
from os.path import isdir, isfile, join

dir_name = './small_asteroid_lightcurve_CFHT_data/'

file_names = [dir_name + f for f in os.listdir(dir_name) if isfile(join(dir_name,f)) ]
# print(file_names)
directory = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

'''
ex fits info

Filename: 1851148p.fits.fz
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU     466   ()      
  1  ccd00         0 CompImageHDU    502   (2112, 4644)   int16   
  2  ccd01         1 CompImageHDU    502   (2112, 4644)   int16   
  3  ccd02         2 CompImageHDU    502   (2112, 4644)   int16   

...and so on...

'''



cts = 0

for f in file_names:
	try:
		file = fits.open ( f )

		# print(f)
	except Exception as e:
		# print('broke at ', f , e)
		continue

	f_id = f.split('/')[-1].split('p')[0]
	
	ind = -1
	for i in range(len(input_file)):
		if f_id in input_file[i,0]:
			ind = i
			break
	dir_ind = -1
	for i in range(len(dir_names)):
		if input_file[ind , 1].split()[1] in dir_names[i]:
			dir_ind = i
			break

	if ind == -1 or dir_ind==-1:
		continue
		
	print( 'file id: ', f , input_file[ind , 1] , input_file[ind , 10] , input_file [ ind , 9 ] )
	print(dir_names[dir_ind])

	chip_id = input_file[ind , 10]

	chips_i_want = []

	# +1 to actual chip numbers for indexing reasons
	if   chip_id=='04' : chips_i_want = [4,6,13,14,15]
	elif chip_id=='13' : chips_i_want = [4,5,6,13,15,22,23,24]
	elif chip_id=='22' : chips_i_want = [13,15,22,24,31,32,33]
	elif chip_id=='31' : chips_i_want = [22,23,24,31,33]

	# fits_header = file[0].header
	primary_hdu = file[0]

	try:
		fits_images = [ file[i].data for i in chips_i_want ]
	except Exception as e:
		print(e , chips_i_want)
		continue
	print(len(fits_images) )

	for i in range(len(chips_i_want)) : 
		c     = chips_i_want[i]
		image = fits_images[i]

		save_to = dir_names[dir_ind] + input_file [ ind , 9 ] + 'on' + str(c) + '.fits'
		print(save_to)

		fits.writeto ( save_to , image , header=file[i+1].header , overwrite=True  )
	# print( i.shape for i in fits_images )


	print()

	
