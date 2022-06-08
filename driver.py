import subprocess, sys
import numpy as np


f = np.loadtxt('star_parameters.csv', delimiter=',', dtype=object, skiprows=1)

L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)

f = f[np.where(L>0)]
L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)
print(len(f))

# print(f)
# index = sys.argv[1]
#for index in range(len(f)-1, -1, -1):
for index in range(21, len(f)):

	print(f[index])
	subprocess.run(['python3', 'magic_star.py', f[index, 0].split(' ')[1], str(L[index]), str(a[index])])#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
