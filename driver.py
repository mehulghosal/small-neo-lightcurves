import subprocess, sys
import numpy as np


f = np.loadtxt('star_parameters.csv', delimiter=',', dtype=object, skiprows=1)

L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)

f = f[np.where(L>0)]
L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)
print(f)

print()

# print(f)
# index = sys.argv[1]
#for index in range(len(f)-1, -1, -1):
for index in range( len(f)):
	#if  ('EV84' in f[index][0]) or ('GE1' in f[index][0]) or ('FF14' in f[index][0]): continue
	#if index == 8: break
	print(f[index])
	print([ 'python3' , 'magic_star.py' , f[index, 0].split(' ')[1] , str(L[index]) , str(a[index]) , str(True) ])
	subprocess.run([ 'python3' , 'magic_star.py' , f[index, 0].split(' ')[1] , str(L[index]) , str(a[index]) , str(True) ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

