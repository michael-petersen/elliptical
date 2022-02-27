

import numpy as np


class SOEllipse(object):
    '''
    #
    # Conic Ellipse fitter
    #     exploiting the quadratic curve nature of the ellipse
    #

    advantages: fast

    disadvantages: does not have flexibility

    notes
    --------
    1. This does, in fact, follow the notation of the wikipedia page
    'Matrix representation of conic sections'

    '''
    @staticmethod
    def fitEllipse(x,y):
        #
        # Take a set of x,y points at fit an ellipse to it
        #
        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
        E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
        n = np.argmax(np.abs(E))
        a = V[:,n]

        return a

    @staticmethod
    def ellipse_center(a):
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num

        return np.array([x0,y0])

    @staticmethod
    def ellipse_angle_of_rotation( a ):
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]

        return 0.5*np.arctan(2*b/(a-c))

    @staticmethod
    def ellipse_axis_length( a ):
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))

        # 27 Feb: I'm putting in abs here. is this bad?
        res1=np.sqrt(np.abs(up/down1))
        res2=np.sqrt(np.abs(up/down2))

        return np.array([res1, res2])




class Ellipse():
    '''
    ellipse definitions for fitting
    '''

    @staticmethod
    def free_ellipse(th,a,b,c):
        '''
        returns generalized ellipse in polar coordinates

        for bar determination following Athanassoula 1990

        '''
        xcomp = ( abs(np.cos(th))**c) / a**c
        ycomp = ( abs(np.sin(th))**c) / b**c
        gell =  ( (xcomp + ycomp) )**(-1./c)
        return gell


    @staticmethod
    def fixed_ellipse(th,a,b):
        '''
        returns c=2 ellipse in polar coordinates


        '''
        xcomp = (( abs(np.cos(th))**2.0) / a**2.0 )
        ycomp = (( abs(np.sin(th))**2.0) / b**2.0 )
        gell =  ( (xcomp + ycomp) )**(-1./2.0)
        return gell


    @staticmethod
    def inside_ellipse(X,Y,A,B,C,rot=0.):
        '''
        inside_ellipse
        determine whether a set of points is inside of an ellipse

        # only tests in first quadrant for power safety



        '''
        rX,rY = X*np.cos(rot)-Y*np.sin(rot),-X*np.sin(rot)-Y*np.cos(rot)

        ellipse_radius = ((abs(rX)/A)**C + (abs(rY)/B)**C)

        yes_ellipse = np.where(ellipse_radius < 1.0)[0]
        ellipse_array = np.zeros(len(X))
        ellipse_array[yes_ellipse] = 1
        return ellipse_array
