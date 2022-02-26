"""
this file shows how to take particle data and make ellipse maps
using exptool (required!)

"""
import numpy as np
import pkg_resources
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# bring in the package itself: check all bugs

from elliptical.trace import map_ellipses

# to build an image from EXP, we'll need exptool
from exptool.io          import particle
from exptool.observables import transform
from exptool.analysis    import pattern
from exptool.observables import visualize


gridsize=128
cres=24
face_extents=0.06
edge_extents=0.02
slice_width=0.1
ktype='gaussian'
npower=5.

indir = ''
runtag = ''

PSPDump = Input(indir+'OUT.'+runtag+'.{0:05d}'.format(300),'star',legacy=True)
PSPDumpt = pattern.BarTransform(PSPDump)

PSPDump = transform.rotate_points(PSPDumpt,0.,0.,0.)
PSPDump.mass = PSPDumpt.mass


X,Y,Z,\
kdeZYz1,kdeZYy1,ZY1,\
kdeXZx1,kdeXZz1,XZ1,\
levels1,levels_edge1 = visualize.kde_pos(PSPDump,\
                                                 gridsize=gridsize,cres=cres,\
                                                 face_extents=face_extents,\
                                                 edge_extents=edge_extents,\
                                                 slice_width=slice_width,\
                                                 ktype=ktype,npower=npower)


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
