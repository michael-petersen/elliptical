"""
different measures for ellipses


references:
Muñoz-Mateos et al. (2013) https://ui.adsabs.harvard.edu/abs/2013ApJ...771...59M/abstract

"""

import numpy as np

class measureEllipse(object):
    '''take a dictionary of ellipses (sorted by semi-major axis length) and return a smattering of measurements.


    '''
    def __init__(self,M,method='pachange'):
        '''constructor. calls out for other measurements.

        M :

        '''

        # recast the parameters as lists
        nlevels = len(M.keys())

        # record the parameters used in measurements
        self.params = dict()

        # prep the arrays
        sma = np.zeros(nlevels)
        ecc = np.zeros(nlevels)
        phi = np.zeros(nlevels)
        indx= np.zeros(nlevels)

        for ik,k in enumerate(M.keys()):
            sma[ik] = M[k]['a']
            ecc[ik] = M[k]['e']
            phi[ik] = M[k]['p']
            indx[ik]= k

        # check that the sorting is correct: a problem if not?
        asort = sma.argsort()

        # ignore the sorting at this point: keep the (inverse) isophotal order
        self.sma = sma[::-1]#[asort]
        self.ecc = ecc[::-1]#[asort]
        self.phi = phi[::-1]#[asort]
        self.indx= indx[::-1]#[asort]

        # now we are ready to compute lengths
        self.ellipdrop      = self._ellip_drop()
        self.maxellip       = self._max_ellip()
        self.localellipmin  = self._first_ellip_min()
        self.ellipdroplimit = self._ellip_drop_below()
        self.maxellipdrop   = self._max_ellip_drop()
        self.ellipchange    = self._ellip_change()
        self.pachange       = self._pa_change()


        if method == 'pachange':
            best_ellipse = np.where(sma==self.pachange)[0][0]
        elif method == 'localellipmin':
            best_ellipse = np.where(sma==self.localellipmin)[0][0]
        elif method == 'ellipchange':
            best_ellipse = np.where(sma==self.ellipchange)[0][0]
        else:
            best_ellipse = np.where(sma==self.maxellip)[0][0]

        # necessary optional logic check
        #print(best_ellipse,self.pachange,M[self.indx[best_ellipse]]['a'],M[best_ellipse]['a'])

        self.bestellipse = M[best_ellipse]


    def print_diagnostics(self):
        '''print all the measurements

        '''
        print('{0:5.4f} | maximum ellipticity'.format(self.maxellip))
        print('{0:5.4f} | position angle change (dphi={1})'.format(self.pachange,self.params['pa_change']))
        #print('{0:5.4f} | sequential ellipticity drop method (de={1})'.format(self.ellipdrop,self.params['ellipdrop']))
        print('{0:5.4f} | ellipticity drop below threshold (e_thresh={1})'.format(self.ellipdroplimit,self.params['ell_threshold']))
        print('{0:5.4f} | maximum ellipticity drop'.format(self.maxellipdrop))
        print('{0:5.4f} | specified ellipticity change (de={1})'.format(self.ellipchange,self.params['el_change']))
        print('{0:5.4f} | local ellipse minimum'.format(self.localellipmin))


    def _ellip_drop(self,drop=0.4):
        '''
        given a list of axis lengths, calculate the length of the bar based on some specified ellipticity drop

        '''

        self.params['ellipdrop'] = drop

        found = False
        j = 2

        while found==False:

            # compute the different in eccentricity from one bin to the next
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

        self.params['ell_threshold'] = drop

        # find all the bins where the ellipticity is above some value (excluded)
        lessthan = np.where( self.ecc >= drop )[0]

        # find the minimum nonzero index
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
        return self.sma[ np.where(np.min(edrop)==edrop)[0][0]]

    def _max_ellip(self):
        """location of maximum ellipticity

        proposed in Munoz-Mateos+ 2013, S3.4 (#1)
        this is treated as the minimum bar length.
        """
        return self.sma[np.nanargmax(self.ecc)]

    def _first_ellip_min(self):
        """location of the first local minimum in the ellipticity after the maximum

        proposed in Munoz-Mateos+ 2013, S3.4 (#2)
        with change=0.1 being the Munoz-Mateos+ 2013 value.


        """

        ellip_index = np.nanargmax(self.ecc)
        max_ellip_value = np.nanmax(self.ecc)

        echange = np.ediff1d(self.ecc,to_end=0.)

        # can put in some tolerance to avoid little wiggles
        while (echange[ellip_index] < 0.):
            ellip_index += 1

        return self.sma[ellip_index-1]


    def _ellip_change(self,change=0.1):
        """location where ellipticity first changes by some amount from the maximum

        proposed in Munoz-Mateos+ 2013, S3.4 (#3)
        with change=0.1 being the Munoz-Mateos+ 2013 value.


        """

        self.params['el_change'] = change

        ellip_index = np.nanargmax(self.ecc)
        max_ellip_value = np.nanmax(self.ecc)

        ellip_diff = 0.

        while (ellip_diff < change):
            ellip_index += 1
            ellip_diff = np.abs(self.ecc[ellip_index] - max_ellip_value)

        return self.sma[ellip_index-1]


    def _pa_change(self,change=10.):
        """location where the angle differs from the position angle at the maximum ellipticity

        proposed in Munoz-Mateos+ 2013, S3.4 (#4)
        with 10 degrees being the Munoz-Mateos+ 2013 value.

        also recently used in JWST bar sample measurements

        inputs
        ------------
        self
        change    : the position angle change tolerance, in degrees.

        """
        # change must be in degrees

        self.params['pa_change'] = change

        ellip_index = np.nanargmax(self.ecc)
        pa_value = self.phi[ellip_index]
        pa_diff = 0.

        while (pa_diff < change):
            ellip_index += 1

            # guard against a possible failure mode
            if ellip_index == len(self.phi):
                print('PA change method failed.')
                return self.sma[ellip_index-1]

            pa_diff = np.abs(self.phi[ellip_index] - pa_value)*180./np.pi

        return self.sma[ellip_index-1]
