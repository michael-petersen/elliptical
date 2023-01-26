"""
trace

definitions to trace the ellipses on an image

follow_contour
make_ellipse
map_ellipses


"""
# standard library
import numpy as np

# the contour-finder. needed!
from skimage.measure import find_contours

# the ellipse definitions
from .ellipse import SOEllipse
from .measure import measureEllipse



def follow_contour(X,Y,Z,level,verbose=0):
    """follow a contour from a given 2d image

    inputs
    ------------
    X          : (2d array) array of X values for image
    Y          : (2d array) array of Y values for image
    Z          : (2d array) surface density values at (X,Y)
    level      : (float)    the contour level to follow
    verbose    : (int)      flag for report

    returns
    ------------
    XCON       : (1d array) x values of the desired contour
    YCON       : (1d array) y values of the desired contour
    """

    # make index boundaries
    uxvals = np.unique(X) # np.unique(X[0,:])
    uxvals = uxvals[np.argsort(uxvals)]
    dx = uxvals[1] - uxvals[0]
    xmin = np.min(uxvals)

    uyvals = np.unique(Y) # np.unique(Y[:,0])
    dy = uyvals[1] - uyvals[0]
    ymin = np.min(uyvals)

    # trace the contour
    res = find_contours(Z,level)

    # extract the contour values
    xcon = []
    ycon = []
    try:
        for i in range(0,len(res[0])):
            xcon.append(dx*res[0][i][0] + xmin)
            ycon.append(dy*res[0][i][1] + ymin)
    except:
        if verbose > 1:
            print('elliptical.trace.follow_contour: No contour found at {}'.format(level))

    # convert the lists to arrays
    XCON = np.array(xcon)
    YCON = np.array(ycon)

    # return the two arrays
    return YCON,XCON



def make_ellipse_conic(xcontours,ycontours):
    """use parametric conic to get basic parameters

    inputs
    -------------
    xcontours     : (1d array) x values of the ellipse points
    ycontours     : (1d array) y values of the ellipse points

    returns
    -------------
    a             : (float) semi-major axis of ellipse
    b             : (float) semi-minor axis of ellipse
    phi           : (float) ellipse angle relative to y=0 axis
    xcenter       : (float) the x centre of the ellipse
    ycenter       : (float) the y centre of the ellipse
    """

    # do the ellipse fit
    ell             = SOEllipse.fitEllipse(xcontours,ycontours)

    # extract parameters
    phi             = SOEllipse.ellipse_angle_of_rotation(ell)
    xcenter,ycenter = SOEllipse.ellipse_center(ell)
    alength         = SOEllipse.ellipse_axis_length(ell)

    # set convention: a is always larger than b
    a = np.max(alength)
    b = np.min(alength)

    # if the second length is larger (e.g. natural b), need to add pi/2 and reverse the centres
    if alength[1] > alength[0]:
        #print('elliptical.make_ellipse_conic: misaligned bar axis')
        phi += np.pi/2.

    # return: cast away any imaginary parts that mistakenly appeared
    return a,b,np.real(phi),np.real(ycenter),np.real(xcenter)


def make_ellipse_parametric(xcontours,ycontours):
    """use parametric ellipse equation to get basic parameters

    inputs
    -------------
    xcontours     : (1d array) x values of the ellipse points
    ycontours     : (1d array) y values of the ellipse points

    returns
    -------------
    a             : (float) semi-major axis of ellipse
    b             : (float) semi-minor axis of ellipse
    phi           : (float) ellipse angle relative to y=0 axis
    xcenter       : (float) the x centre of the ellipse
    ycenter       : (float) the y centre of the ellipse
    """

    # do the ellipse fit
    ell             = SOEllipse.fitEllipse(xcontours,ycontours)

    # extract parameters
    phi             = SOEllipse.ellipse_angle_of_rotation(ell)
    xcenter,ycenter = SOEllipse.ellipse_center(ell)
    alength         = SOEllipse.ellipse_axis_length(ell)

    # set convention: a is always larger than b
    #   this doesn't affect the angles, right?
    a = np.max(alength)
    b = np.min(alength)

    # return: cast away any imaginary parts that mistakenly appeared
    return a,b,np.real(phi),np.real(xcenter),np.real(ycenter)




def map_ellipses(X,Y,Z,minZ,maxZ,numZ=16,CENTERTOL=1.,PHITOL=7.,ETOL=0.5,optimal=False,verbose=0,method='pachange'):
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
    CENTERTOL  : (float)    tolerance distance an ellipse may range from the centre
    PHITOL     : (float)    tolerance (radian) angle for defining ellipses (used if bar is pre-aligned)
    ETOL       : (float)    tolerance for elliptical-ness in defining best ellipse
    optimal    : (bool)     if True, return only the best-fit ellipse
    verbose    : (int)      verbosity flag. Increase for more report.
    method: (string)   method to use for designating the best-fit ellipse

    returns
    -----------
    M          : (dict)
      if optimal, returns a one-level dictionary with x,y of the best-fit ellipse, and a (the semi-major axis)
      if !optimal, returns a two-level dictionary with all drawn ellipses.

    """

    # allocate a storage location
    M = dict()

    # set up angular samples
    th = np.arange(0,2*np.pi, 0.01)

    # define the contour levels to try drawing
    ctestvals = np.linspace(minZ,maxZ,numZ)

    # initialise variables
    maxa      = 0.
    best_cval = 0.
    cnum      = 0

    previousa = 1.e6

    # loop through contour levels
    for cval in ctestvals:

        try:
            # trace the contour
            XCON,YCON = follow_contour(X,Y,Z,cval,verbose=verbose)

            # make the ellipse from the countour
            #a,b,phi,xcenter,ycenter = make_ellipse_conic(XCON,YCON)

            # contours come out in reverse order
            a,b,phi,xcenter,ycenter = make_ellipse_conic(YCON,XCON)

            # if a good ellipse, save values
            if ((np.sqrt(xcenter*xcenter + ycenter*ycenter) < CENTERTOL) & (phi < PHITOL)):

                M[cnum] = dict()
                xx = xcenter + a*np.cos(th)*np.cos(phi) - b*np.sin(th)*np.sin(phi)
                yy = ycenter + a*np.cos(th)*np.sin(phi) + b*np.sin(th)*np.cos(phi)

                # save ellipse parameters
                M[cnum]['x'] = xx
                M[cnum]['y'] = yy
                M[cnum]['a'] = a
                M[cnum]['b'] = b
                M[cnum]['e'] = 1.-b/a
                M[cnum]['p'] = phi
                M[cnum]['l'] = cval
                M[cnum]['xc'] = xcenter
                M[cnum]['yc'] = ycenter

                # advance the ellipse counter
                cnum += 1

        # if the ellipse drawing fails, move on
        except:
            pass

    if optimal:
        ME = measureEllipse(M,method=method)

        return ME.bestellipse

    else:
        if (verbose>0):
            print("You requested {} levels, found {} valid levels.".format(numZ,cnum-1))
        return M
