# different measures for ellipses

import numpy as np

class measureEllipse(object):
    '''take a dictionary of ellipses (sorted by semi-major axis length) and return a smattering of measurements.


    '''
    def __init__(self,M):
        '''constructor. calls out for other measurements.'''

        # recast the parameters as lists
        nlevels = len(M.keys())

        sma = np.zeros(nlevels)
        ecc = np.zeros(nlevels)
        phi = np.zeros(nlevels)

        for ik,k in enumerate(M.keys()):
            sma[ik] = M[k]['a']
            ecc[ik] = M[k]['e']
            phi[ik] = M[k]['p']

        # check that the sorting is correct
        asort = sma.argsort()
        self.sma = sma[asort]
        self.ecc = ecc[asort]
        self.phi = phi[asort]

        # now we are ready to compute lengths
        self.ellipdrop      = self._ellip_drop()
        self.maxellip       = self._max_ellip()
        self.ellipdroplimit = self._ellip_drop_below()
        self.maxellipdrop   = self._max_ellip_drop()
        self.ellipchange    = self._ellip_change()
        self.pachange       = self._pa_change()

    def print_diagnostics(self):
        '''print all the measurements'''
        print('{0:5.4f} | maximum ellipticity'.format(self.maxellip))
        print('{0:5.4f} | position angle change'.format(self.pachange))
        print('{0:5.4f} | sequential ellipticity drop method'.format(self.ellipdrop))
        print('{0:5.4f} | ellipticity drop below threshold'.format(self.ellipdroplimit))
        print('{0:5.4f} | maximum ellipticity drop'.format(self.maxellipdrop))
        print('{0:5.4f} | specified ellipticity change'.format(self.ellipchange))


    def _ellip_drop(self,drop=0.4):
        '''
        given a list of axis lengths, calculate the length of the bar based on some specified ellipticity drop

        '''
        found = False
        j = 2
        while found==False:
            d = self.ecc[j-1] - self.ecc[j]
            if d > drop:
                found = True
                print('INDEX VALUE is {0:d}'.format(j-1))
            j += 1
            if j==len(self.sma):
                found = True
                j=2
        return self.sma[j-2]

    def _ellip_drop_below(self,drop=0.4):
        '''where does the ellipticity first drop below some value?

        '''
        lessthan = np.where( self.ecc >= drop )[0]
        if len(lessthan) > 0:
            if np.max(lessthan) < 1:
                minbin = 1
            else:
                minbin = np.max(lessthan)
        else:
            minbin = 1
        return self.sma[minbin]

    def _max_ellip_drop(self):
        """location of largest eccentricity change"""
        edrop = np.ediff1d(self.ecc,to_end=0.)
        return self.sma[ np.where(np.min(edrop)==edrop)[0]]

    def _max_ellip(self):
        """location of maximum ellipticity

        proposed in Munoz-Mateos+ 2013, S3.4
        this is treated as the minimum bar length.
        """
        return self.sma[np.nanargmax(self.ecc)]

    def _ellip_change(self,change=0.1):
        """location where ellipticity first changes by some amount

        proposed in Munoz-Mateos+ 2013, S3.4
        with change=0.1 being the Munoz-Mateos+ 2013 value.
        """
        ellip_index = np.nanargmax(self.ecc)
        max_ellip_value = np.nanmax(self.ecc)
        ellip_diff = 0.
        while (ellip_diff < change):
            ellip_index += 1
            ellip_diff = np.abs(self.ecc[ellip_index] - max_ellip_value)
        return self.sma[ellip_index-1]


    def _pa_change(self,change=10.):
        """location where the angle differs from the position angle at the maximum ellipticity

        proposed in Munoz-Mateos+ 2013, S3.4
        with 10 degrees being the Munoz-Mateos+ 2013 value.

        inputs
        ------------
        self
        change    : the position angle change tolerance, in degrees.

        """
        # change must be in degrees

        ellip_index = np.nanargmax(self.ecc)
        pa_value = self.phi[ellip_index]
        pa_diff = 0.
        while (pa_diff < change):
            ellip_index += 1
            pa_diff = np.abs(self.phi[ellip_index] - pa_value)*180./np.pi
        return self.sma[ellip_index-1]
