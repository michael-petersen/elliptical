"""
this file tests the basic ellipse-drawing capabilities

"""
import numpy as np
import pkg_resources
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# bring in the package itself: check all bugs
import elliptical
print('elliptical version {}'.format(elliptical.__version__))

from elliptical.trace import map_ellipses

from elliptical.measure import measureEllipse

from elliptical.ellipse import inside_ellipse

# identify the testing files
g1 = pkg_resources.resource_filename('elliptical','data/galaxy2.dat')

# unpack a test image
nmodel1 = np.genfromtxt(g1,max_rows=1)
xdim,ydim = int(nmodel1[0]),int(nmodel1[1])

model1 = np.genfromtxt(g1,skip_header=1)
X,Y,Z = model1[:,0].reshape([xdim,ydim]),model1[:,1].reshape([xdim,ydim]),model1[:,2].reshape([xdim,ydim])


plt.figure(facecolor='white')

# draw the image
plt.contourf(X,Y,Z,24,cmap=cm.inferno)

# map ellipses (all)
M = map_ellipses(X,Y,Z,-6.5,-4.,numZ=64,verbose=1)

for k in M.keys():
    plt.plot(M[k]['x'],M[k]['y'],color='black',lw=1.)

    # prove the de-rotation works
    #xell = M[k]['x']*np.cos(M[k]['p']) + M[k]['y']*np.sin(M[k]['p'])
    #yell =-M[k]['x']*np.sin(M[k]['p']) + M[k]['y']*np.cos(M[k]['p'])
    #plt.plot(xell,yell,color='white',lw=1.)



# map ellipse (best-fit only)
MB = map_ellipses(X,Y,Z,-6.5,-4.,numZ=64,optimal=True,method='pachange')
plt.plot(MB['x'],MB['y'],color='white',linestyle='dashed',lw=1.)

print(MB['p'])
print(MB['a'],MB['b'])

mask = inside_ellipse(MB['a'],MB['b'],MB['p'],MB['xc'],MB['yc'],X,Y)
#plt.contourf(X,Y,Z*mask,24,cmap=cm.inferno)


plt.xlabel('X [scale lengths]')
plt.ylabel('Y [scale lengths]')
plt.tight_layout()
plt.savefig('galaxymodel1.png')


ME = measureEllipse(M)
#print(ME.maxellip)

plt.figure(facecolor='white')
plt.plot(ME.sma,ME.phi,color='black',drawstyle='steps-mid')
plt.plot([ME.maxellip,ME.maxellip],[0.,1.],color='grey',linestyle='dashed',lw=1.)
#plt.plot([MB['a'],MB['a']],[0.,1.],color='red',linestyle='dashed',lw=1.)

plt.axis([0.,np.nanmax(X),0.,1.])
plt.axis([0.,np.nanmax(X),-2,2.])

plt.xlabel('a [scale lengths]')
plt.ylabel('ecc')
plt.tight_layout()
plt.savefig('galaxymodel1measure.png')
