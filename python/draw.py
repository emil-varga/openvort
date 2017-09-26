import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

points0 = np.loadtxt('../data/step0000.dat')
points0a = np.loadtxt('../data/step0250.dat')
points1 = np.loadtxt('../data/step0500.dat')

pointsa = np.loadtxt('../v1.dat')
pointsb = np.loadtxt('../v2.dat')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.plot(pointsa[:,1], pointsa[:,2], pointsa[:,3])
#ax.plot(pointsb[:,1], pointsb[:,2], pointsb[:,3])
#ax = fig.add_subplot(111);

#ax.plot(points0[:,1], points0[:,3], '-o', markersize=10)
#ax.plot(pointsa[:,1], pointsa[:,3], '-ro', markersize=5)
#ax.set_aspect('equal')
ax.plot(points0[:,1], points0[:,2], points0[:,3])
ax.plot(points0a[:,1], points0a[:,2], points0a[:,3])
ax.plot(points1[:,1], points1[:,2], points1[:,3])

f, ax = plt.subplots(1,1)
ax.plot(points0[:,1], points0[:,2], '-o')
ax.plot(points0a[:,1], points0a[:,2], '-o')
ax.plot(points1[:,1], points1[:,2], '-o')

ers = []
ics = range(500)
for n in ics:
    file = "../data/step{:04d}.dat".format(n)
    d = np.loadtxt(file)

    error = (d[:,1] - points0[:,1])**2 + (d[:,2] - points0[:,2])**2
    ers.append(error.sum())

f, ax = plt.subplots(1,1)
ax.semilogy(ics, ers)
plt.show()
