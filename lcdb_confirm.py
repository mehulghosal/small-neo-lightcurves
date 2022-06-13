import subprocess, sys
import numpy as np


f = np.loadtxt('star_parameters.csv', delimiter=',', dtype=object, skiprows=1)

L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)

# f = f[np.where(L>0)]
# L = np.array(f[:,1], dtype=float)
# a = np.array(f[:,2], dtype=float)
# print(len(f))

lcdb_summary = np.loadtxt('/home/mehul/code/ast-lightcurve-database_V4_0/data/lc_summary.csv', delimiter=',', skiprows=22, dtype=object)

# print(f)
# index = sys.argv[1]
#for index in range(len(f)-1, -1, -1):
for index in range(len(f)):

	# print(f[index])
	obj_name = f[index, 0].split(' ')
	# print(obj_name[0].join(obj_name[1]))
	obj_name = obj_name[0] + " " + obj_name[1]
	# print(obj_name)
	lcdb_ind = -1
	for j in range(len(lcdb_summary)):

		if obj_name in lcdb_summary[j, 1]: 
			lcdb_ind = j
			break
	# print(lcdb_ind)
	if not lcdb_ind == -1: 

		print(obj_name, lcdb_summary[lcdb_ind, 19], lcdb_summary[lcdb_ind, 8])
		# print(lcdb_summary[lcdb_ind])
		# pass

	# subprocess.run(['python3', 'magic_star.py', f[index, 0].split(' ')[1], str(L[index]), str(a[index])])#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
