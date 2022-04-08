import subprocess, sys
import numpy as np


f = np.loadtxt('star_parameters.csv', delimiter=',', dtype=object, skiprows=1)

L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)

f = f[np.where(L>0)]
L = np.array(f[:,1], dtype=float)
a = np.array(f[:,2], dtype=float)
print(len(f))


index = sys.argv[1]
subprocess.run(['python3', 'magic_star.py', f[index, 0].split(' ')[1], L[index], a[index]])