import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.close('all')

#vi, pos, vel, s', s''

points0 = np.loadtxt('../data_rec/step0000.dat')
points1 = np.loadtxt('../data_rec/step0100.dat')
points2 = np.loadtxt('../data_rec/step0200.dat')

f, ax = plt.subplots(1,1)
ax.semilogy(np.abs(points1[:,6]-points0[:,6]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(points0[:,1], points0[:,2], points0[:,3])
ax.plot(points1[:,1], points1[:,2], points1[:,3])
ax.plot(points2[:,1], points2[:,2], points2[:,3])

f, ax = plt.subplots(1,1)
ax.plot(points0[:,1], points0[:,3], '-o')
ax.plot(points1[:,1], points1[:,3], '-o')
ax.plot(points2[:,1], points2[:,3], '-o')

ds = []

for n in range(points0.shape[0] - 1):
    p0 = points0[n,10:]
    p1 = points0[n+1, 4:7]
    
    ds.append(((p0)**2).sum())
f, ax = plt.subplots(1,1)
ax.plot(ds, '-o')

#ers = []
#ics = range(500)
#for n in ics:
#    file = "../data_rk/step{:04d}.dat".format(n)
#    d = np.loadtxt(file)
#
##    error = (d[:,1] - points0[:,1])**2 + (d[:,2] - points0[:,2])**2
#    error = (d[:,3] - d[:,3].mean())**2
#    ers.append(error.sum())
#
#f, ax = plt.subplots(1,1)
#ax.plot(ics, ers)
#plt.show()
