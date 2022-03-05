"""deproject

theoretical considerations for deprojecting ellipses

see Appendix A of Gadotti et al. (2007) for derivations
https://ui.adsabs.harvard.edu/abs/2007MNRAS.381..943G/abstract

"""

import numpy as np

def st_from_xy(s,t,x,y,alpha):
    return x*np.cos(alpha) + y*np.sin(alpha),y*np.cos(alpha)-x*np.sin(alpha)

class Deproject(object):

    def __init__(self,a,b,alpha,i):

        self.a     = a
        self.b     = b
        self.alpha = alpha
        self.i     = i

        A = self._A()
        B = self._B()
        C = self._C()
        D = self._D()
        F = self._F()
        G = self._G()

        s1 = Deproject._s1(A,B,C,D,F,G)
        s2 = Deproject._s2(A,B,C,D,F,G)

        self.sma = s1#np.nanmax([s1,s2])
        self.smb = s2#np.nanmin([s1,s2])
        self.ecc = 1. - (self.smb/self.sma)

    def _A(self):
        return np.cos(self.alpha)**2./(self.a**2.) + np.sin(self.alpha)**2./(self.b**2.)

    def _B(self):
        cossinalpha = np.cos(self.alpha)*np.sin(self.alpha)
        Bunprime = (cossinalpha/(self.a**2.) - cossinalpha/(self.b**2.))
        return Bunprime*np.cos(self.i)

    def _C(self):
        Cunprime = np.sin(self.alpha)**2./(self.a**2.) + np.cos(self.alpha)**2./(self.b**2.)
        return Cunprime*np.cos(self.i)**2.

    def _D(self):
        return 0.

    def _F(self):
        return 0.

    def _G(self):
        return -1.

    @staticmethod
    def _s1(A,B,C,D=0.,F=0.,G=-1.):
        numerator = 2*(A*F*F + C*D*D + G*B*B - 2*B*D*F - A*C*G)
        denominator = (B*B-A*C)*((C-A)*np.sqrt(1+4*B*B/((A-C)**2.))-(C+A))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s1_simple(A,B,C):
        numerator = 2*((-1)*B*B + A*C)
        denominator = (B*B-A*C)*((C-A)*np.sqrt(1+4*B*B/((A-C)**2.))-(C+A))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s2(A,B,C,D=0.,F=0.,G=-1.):
        numerator = 2*(A*F*F + C*D*D + G*B*B - 2*B*D*F - A*C*G)
        denominator = (B*B-A*C)*((A-C)*np.sqrt(1+4*B*B/((A-C)**2.))-(C+A))
        return np.sqrt(numerator/denominator)

    @staticmethod
    def _s2_simple(A,B,C):
        """using the known D,F,G values)"""
        numerator = 2*((-1)*B*B + A*C)
        denominator = (B*B-A*C)*((A-C)*np.sqrt(1+4*B*B/((A-C)**2.))-(C+A))
        return np.sqrt(numerator/denominator)
