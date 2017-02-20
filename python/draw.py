import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

points0 = np.loadtxt('step0.dat')
points1 = np.loadtxt('step1.dat')
pointsa = np.loadtxt('steapa.dat')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax = fig.add_subplot(111);

#ax.plot(points0[:,1], points0[:,3], '-o', markersize=10)
#ax.plot(pointsa[:,1], pointsa[:,3], '-ro', markersize=5)
#ax.set_aspect('equal')
ax.plot(points0[:,1], points0[:,2], points0[:,3])
ax.plot(points1[:,1], points1[:,2], points1[:,3])

plt.show()
