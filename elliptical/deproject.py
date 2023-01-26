"""deproject

theoretical considerations for deprojecting ellipses

see Appendix A of Gadotti et al. (2007) for derivations
https://ui.adsabs.harvard.edu/abs/2007MNRAS.381..943G/abstract

see additional testing for verification in Zou et al. (2014)
https://ui.adsabs.harvard.edu/abs/2014ApJ...791...11Z/abstract

"""

import numpy as np

def st_from_xy(s,t,x,y,alpha):
    return x*np.cos(alpha) + y*np.sin(alpha),y*np.cos(alpha)-x*np.sin(alpha)

class Deproject(object):

    def __init__(self,a,b,alpha,i):
        """

        This is presumably solvable in Mathematica. It would be nice to check the results.

        the ellipse needs to be pre-centred.

        We have some ellipse, traced at points x and y.
        Ideally we also know the rotation angle (the angle between the bar and the line of nodes) alpha.
        And the inclination angle i.

        **Note that equations for s1 and s2 have singularities when i = 0 and α = ±nπ/4 (n being a positive integer). A singularity also appears in equation the position angle of the deprojected ellipse when α = 0, ±nπ/2. If i = π/2, the above three equations diverge.

        """

        self.a     = a
        self.b     = b
        self.alpha = alpha
        self.i     = i

        A = self._Aprime()
        B = self._Bprime()
        C = self._Cprime()
        D = self._Dprime()
        F = self._Fprime()
        G = self._Gprime()

        s1 = Deproject._s1(A,B,C,D,F,G)
        s2 = Deproject._s2(A,B,C,D,F,G)

        self.sma = np.nanmax([s1,s2])
        self.smb = np.nanmin([s1,s2])
        self.ecc = 1. - (self.smb/self.sma)

        # the ellipse comes out by defininition with phi=0, xcentre=ycentre=0

    def _A(self):
        """equation A4"""
        return (np.cos(self.alpha)**2.)/(self.a**2.) + (np.sin(self.alpha)**2.)/(self.b**2.)

    def _B(self):
        """equation A5"""
        cossinalpha = np.cos(self.alpha)*np.sin(self.alpha)
        return (cossinalpha/(self.a**2.) - cossinalpha/(self.b**2.)) # this is equation A5

    def _C(self):
        """equation A6"""
        return np.sin(self.alpha)**2./(self.a**2.) + np.cos(self.alpha)**2./(self.b**2.)

    def _D(self):
        """after equation A6"""
        return 0.

    def _F(self):
        """after equation A6"""
        return 0.

    def _G(self):
        """after equation A6"""
        return -1.

    def _Bprime(self):
        """equation A8"""
        return self._B()*np.cos(self.i)

    def _Cprime(self):
        """equation A9"""
        return self._C()*(np.cos(self.i)**2.)

    def _Aprime(self):
        """after equation A9"""
        return self._A()

    def _Dprime(self):
        """after equation A9"""
        return 0.

    def _Fprime(self):
        """after equation A9"""
        return 0.

    def _Gprime(self):
        """after equation A9"""
        return self._G()

    @staticmethod
    def _s1(Aprime,Bprime,Cprime,Dprime=0.,Fprime=0.,Gprime=-1.):
        """equation A10"""
        numerator = 2*(Aprime*Fprime*Fprime + Cprime*Dprime*Dprime + Gprime*Bprime*Bprime - 2*Bprime*Dprime*Fprime - Aprime*Cprime*Gprime)
        denominator = (Bprime*Bprime-Aprime*Cprime)*((Cprime-Aprime)*np.sqrt(1.+4*Bprime*Bprime/((Aprime-Cprime)**2.))-(Cprime+Aprime))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s1_simple(Aprime,Bprime,Cprime):
        """equation A10 with defaults applied"""
        numerator = 2*((-1.)*Bprime*Bprime + Aprime*Cprime)
        denominator = (Bprime*Bprime-Aprime*Cprime)*((Cprime-Aprime)*np.sqrt(1.+4*Bprime*Bprime/((Aprime-Cprime)**2.))-(Cprime+Aprime))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s2(Aprime,Bprime,Cprime,Dprime=0.,Fprime=0.,Gprime=-1.):
        """equation A11"""
        numerator = 2*(Aprime*Fprime*Fprime + Cprime*Dprime*Dprime + Gprime*Bprime*Bprime - 2*Bprime*Dprime*Fprime - Aprime*Cprime*Gprime)
        denominator = (Bprime*Bprime-Aprime*Cprime)*((Aprime-Cprime)*np.sqrt(1.+4*Bprime*Bprime/((Aprime-Cprime)**2.))-(Cprime+Aprime))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s2_simple(Aprime,Bprime,Cprime):
        """using the known D,F,G values)"""
        numerator = 2*((-1.)*Bprime*Bprime + Aprime*Cprime)
        denominator = (Bprime*Bprime-Aprime*Cprime)*((Aprime-Cprime)*np.sqrt(1.+4*Bprime*Bprime/((Aprime-Cprime)**2.))-(Cprime+Aprime))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _deprojected_position_angle(Aprime,Bprime,Cprime):
        """equation A12"""
        return -0.5 * (np.arctan((Cprime-Aprime)/(2*Bprime))**(-1))
