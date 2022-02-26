"""
trace

collection of 



"""

import numpy as np

# the contour-finder. needed!
from skimage.measure import find_contours

# the ellipse definitions
from .ellipse import SOEllipse



def follow_contour(xx,yy,arr,level,verbose=0):
    """follow a contour from a given 2d image

    xx and yy are assumed to be evenly spaced

    """

    dx = np.unique(xx[0,:])[1] - np.unique(xx[0,:])[0]
    xmin = np.min(np.unique(xx[0,:]))

    dy = np.unique(yy[:,0])[1] - np.unique(yy[:,0])[0]
    ymin = np.min(np.unique(yy[:,0]))

    res = find_contours(arr,level)

    xcon = []
    ycon = []
    try:
        for i in range(0,len(res[0])):
            xcon.append(dx*res[0][i][0] + xmin)
            ycon.append(dy*res[0][i][1] + ymin)
    except:
        if verbose:
            print('draw.follow_contour: No countour found at {}'.format(level))

    XCON = np.array(xcon)
    YCON = np.array(ycon)
    return YCON,XCON



def make_ellipse(xcontours,ycontours):
    #
    # use parametric conic to get basic parameters
    #
    ell = SOEllipse.fitEllipse(xcontours,ycontours)
    phi = SOEllipse.ellipse_angle_of_rotation(ell)
    xcenter,ycenter = SOEllipse.ellipse_center(ell)
    alength = SOEllipse.ellipse_axis_length(ell)
    a = np.max(alength)
    b = np.min(alength)
    return a,b,phi,xcenter,ycenter



def map_ellipses(X,Y,Z,minZ,maxZ,numZ=16,centertol=1.,optimal=False):
    """
    create a map of ellipses from an image

    inputs
    -----------
    X          : (2d array) array of X values for image
    Y          : (2d array) array of Y values for image
    Z          : (2d array) surface density values at (X,Y)
    minZ       : (float)    minimum contour level to try drawing
    maxZ       : (float)    maximum contour level to try drawing
    numZ       : (int)      number of ellipses to try and draw
    centertol  : (float)    maximum distance an ellipse may range from the centre
    optimal    : (bool)     if True, return only the best-fit ellipse

    returns
    -----------
    M          : (dict)
      if optimal, returns a one-level dictionary with x,y of the best-fit ellipse, and a (the semi-majro axis)
      if !optimal, returns a two-level dictionary with all drawn ellipses.

    """

    # allocate a storage location
    M = dict()

    R = np.arange(0,2*np.pi, 0.01)

    maxa=0.
    best_cval=0.
    ctestvals = np.linspace(minZ,maxZ,16)
    cnum = 0
    for cval in ctestvals:
        #print(cval)
        try:
            XCON,YCON = follow_contour(X,Y,Z,cval)
            #print(XCON,YCON)
            a,b,phi,xcenter,ycenter = make_ellipse(XCON,YCON)
            #print((xcenter*xcenter + ycenter*ycenter)**0.5)
            if ((xcenter*xcenter + ycenter*ycenter)**0.5 < centertol) & (phi < 0.3):
                M[cnum] = dict()
                xx = xcenter + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
                yy = ycenter + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)
                M[cnum]['x'] = xx
                M[cnum]['y'] = yy
                M[cnum]['a'] = a
                M[cnum]['e'] = 1.-b/a
                M[cnum]['p'] = phi

                cnum += 1

            if (((b/a)<0.5) & (a < np.nanmax(X)) & ((xcenter*xcenter + ycenter*ycenter)**0.5 < centertol)):
                if a>maxa:
                    maxa = a
                    best_cval = cnum - 1
        except:
            pass

    if optimal:
        return M[best_cval]
    else:
        return M
