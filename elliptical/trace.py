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
    dx = np.unique(X[0,:])[1] - np.unique(X[0,:])[0]
    xmin = np.min(np.unique(X[0,:]))

    dy = np.unique(Y[:,0])[1] - np.unique(Y[:,0])[0]
    ymin = np.min(np.unique(Y[:,0]))

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
        if verbose > 0:
            print('elliptical.trace.follow_contour: No contour found at {}'.format(level))

    # convert the lists to arrays
    XCON = np.array(xcon)
    YCON = np.array(ycon)

    # return the two arrays
    return YCON,XCON



def make_ellipse(xcontours,ycontours):
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
    #   this doesn't affect the angles, right?
    a = np.max(alength)
    b = np.min(alength)

    # return
    return a,b,phi,xcenter,ycenter



def map_ellipses(X,Y,Z,minZ,maxZ,numZ=16,CENTERTOL=1.,PHITOL=0.3,ETOL=0.5,optimal=False):
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

    returns
    -----------
    M          : (dict)
      if optimal, returns a one-level dictionary with x,y of the best-fit ellipse, and a (the semi-majro axis)
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

    # loop through contour levels
    for cval in ctestvals:

        try:
            # trace the contour
            XCON,YCON = follow_contour(X,Y,Z,cval)

            # make the ellipse from the countour
            a,b,phi,xcenter,ycenter = make_ellipse(XCON,YCON)

            # if a good ellipse, save values
            if ((np.sqrt(xcenter*xcenter + ycenter*ycenter) < CENTERTOL) &
            (phi < PHITOL)):

                M[cnum] = dict()
                xx = xcenter + a*np.cos(th)*np.cos(phi) - b*np.sin(th)*np.sin(phi)
                yy = ycenter + a*np.cos(th)*np.sin(phi) + b*np.sin(th)*np.cos(phi)

                # save ellipse parameters
                M[cnum]['x'] = xx
                M[cnum]['y'] = yy
                M[cnum]['a'] = a
                M[cnum]['e'] = 1.-b/a
                M[cnum]['p'] = phi

                # advance the ellipse counter
                cnum += 1

            # identify the largest ellipse satisfying bar constraints
            if ((a < np.nanmax(X)) &
            ((b/a)<ETOL) &
            (np.sqrt(xcenter*xcenter + ycenter*ycenter) < CENTERTOL) &
            (phi < PHITOL)):

                # check if this ellipse is larger than the current max
                if a>maxa:
                    maxa = a
                    best_cval = cnum - 1

        # if the ellipse drawing fails, move on
        except:
            pass

    if optimal:
        return M[best_cval]
    else:
        return M
