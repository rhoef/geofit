# -*- coding: utf-8 -*-
"""
ellipse_geo.py

Perform the best fit to a 2d data set i.e.
minimizing the geometrical distance.
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'
__svn_id__ = '$Id$'

from ellipse import EllipseBase
import numpy as np
from scipy.optimize import leastsq


class EllipseGeometric(EllipseBase):

    def __init__(self, x, y, p0=None):
        self.fit()
        self.normal_from()
        self.params()

    def mfunc(self, p, x, y):
        phi = self.phi() - p[4]
        # cx, cy, a, b, alpha, phi
        xbar = p[2]*np.cos(phi)
        ybar = p[3]*np.sin(phi)
#        rotm = self.rotmat(p[4])
        rotx = np.cos(p[4])*xbar - np.sin(p[4])*ybar
        roty = np.sin(p[4])*xbar + np.cos(p[4])*ybar

        print x.shape
        print rotx.shape
        print y.shape
        print roty.shape

#        import pdb; pdb.set_trace()

       # di = np.sqrt((p[0]-x)**2 + (p[1]-y)**2) - np.sqrt(rotx**2 + roty**2)
        di = (x - p[0])**2 + (y -p[1])**2  - (rotx**2 + roty**2)
        return di

    def func(self, p, x, y):
        return (x-p[0])**2/p[2]**2 + (y-p[1])**2/p[3]**2 - 1

    def chi2(self, x, y):
        return lambda p: self.func(p, x, y).flatten()

    def xbar(self, a, b, phi):
        return a*np.cos(phi)

    def rotmat(self, alpha):
        return np.matrix([np.cos(alpha), -np.sin(alpha)],
                         [np.sin(alpha), np.cos(alpha)])

    def fit(self):

        efunc = self.chi2(self.x, self.y)
        p0 = self.p0
        popt, err = leastsq(efunc, p0, full_output=False)
#        import pdb; pdb.set_trace()

        self.center = popt[:2]
        self.a = popt[2]
        self.b = popt[3]
        self.phi0 = np.degrees(popt[4])

    @property
    def p0(self):

        if not self._p0:
        #        p0 = np.append(np.array([0.0, 0.0, 111, 112, 0.0]), phi0)
            return np.append( [0.0, 0.0, 111, 112, 0.0], self.phi)
        else:
            return self._p0
#           return np.append(self._p0, self.phi())

    def phi(self, x0=0.0, y0=0.0):

        return np.arctan2(self.y-y0, self.x-x0)
