import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

points0 = np.loadtxt('step0.dat')
points1 = np.loadtxt('step1.dat')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(points0[:,1], points0[:,2], points0[:,3])
ax.plot(points1[:,1], points1[:,2], points1[:,3])

plt.show()
