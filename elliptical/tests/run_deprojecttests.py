"""
this file tests some deprojection of ellipses

"""
import numpy as np
import pkg_resources
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# bring in the package itself: check all bugs
import elliptical

from elliptical.deproject import Deproject
from elliptical.ellipse import SOEllipse
from elliptical.trace import make_ellipse

D = Deproject(2.,1.,0.,np.pi/4.)
print(D.sma,D.smb,D.ecc)


def construct_tait_bryan(xrotation,yrotation,zrotation):
    radfac = np.pi/180.

    # set rotation in radians
    a = xrotation*radfac # xrotation (the tip into/out of page)
    b = yrotation*radfac # yrotation
    c = zrotation*radfac # zrotation

    # construct the rotation matrix TAIT-BRYAN method
    # (x-y-z, extrinsic rotations)
    Rx = np.array([[1.,0.,0.],[0.,np.cos(a),np.sin(a)],[0.,-np.sin(a),np.cos(a)]])
    Ry = np.array([[np.cos(b),0.,-np.sin(b)],[0.,1.,0.],[np.sin(b),0.,np.cos(b)]])
    Rz = np.array([[np.cos(c),np.sin(c),0.,],[-np.sin(c),np.cos(c),0.],[0.,0.,1.]])
    Rmatrix = np.dot(Rx,np.dot(Ry,Rz))
    return Rmatrix

def rotate_xy(xx,yy,Rmatrix):
    pts = np.array([xx,yy,xx*0.])
    tmp = np.dot(pts.T,Rmatrix)
    xp = tmp[:,0]
    yp = tmp[:,1]
    zp = tmp[:,2]
    return xp,yp

# draw a circle (first test into/out of page)

th = np.linspace(0.,2.*np.pi,200)
rc = 1.
xx = 2.*rc*np.cos(th)
yy = rc*np.sin(th)


plt.figure()

plt.plot(xx,yy,color='black')

yrot = 50.
zrot = 50.

for indx in np.linspace(10.,85.,10):

    # the last rotation has no bearing on the deprojection (or the measured ellipse)
    # the y rotation will need some work, though.
    xp,yp = rotate_xy(xx,yy,construct_tait_bryan(indx,yrot,zrot))
    EE = make_ellipse(xp,yp)


    D = Deproject(EE[0],EE[1],yrot*np.pi/180.,indx*np.pi/180.)
    print(np.round(indx,1),np.round(EE[0],2),np.round(EE[1],2),np.round(D.sma,2),np.round(D.smb,2))

    plt.plot(xp,yp,color=cm.viridis(indx/90.,1.))

plt.xlabel('X [scale lengths]')
plt.ylabel('Y [scale lengths]')
plt.tight_layout()
plt.savefig('deproject1.png')
