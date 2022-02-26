"""
this file tests the basic ellipse-drawing capabilities

"""
import numpy as np
import pkg_resources
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# bring in the package itself: check all bugs
import elliptical

from elliptical.trace import map_ellipses

# identify the testing files
g1 = pkg_resources.resource_filename('elliptical','data/galaxy1.dat')

# unpack a test image
nmodel1 = np.genfromtxt(g1,max_rows=1)
xdim,ydim = int(nmodel1[0]),int(nmodel1[1])

model1 = np.genfromtxt(g1,skip_header=1)
X,Y,Z = model1[:,0].reshape([xdim,ydim]),model1[:,1].reshape([xdim,ydim]),model1[:,2].reshape([xdim,ydim])


plt.figure()

# draw the image
plt.contourf(X,Y,Z,24,cmap=cm.inferno)

# map ellipses (all)
M = map_ellipses(X,Y,Z,-6.5,-4.,numZ=16)

for k in M.keys():
    plt.plot(M[k]['x'],M[k]['y'],color='black',lw=1.)

# map ellipse (best-fit only)
M = map_ellipses(X,Y,Z,-6.5,-4.,numZ=16,optimal=True)
plt.plot(M['x'],M['y'],color='white',linestyle='dashed',lw=1.)

plt.xlabel('X [scale lengths]')
plt.ylabel('Y [scale lengths]')
plt.tight_layout()
plt.savefig('galaxymodel1.png')
