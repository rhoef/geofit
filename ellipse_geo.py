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

def flatzip(x, y):
    zipped = np.array(zip(x, y))
    return zipped.flatten()

class EllipseGaussNewton(EllipseBase):

    def __init__(self, x, y, p0):
        super(EllipseGaussNewton, self).__init__(x, y)
        self.properties = {"circ_tol": 1e-5,
                           "max_iterations": 666,
                           "converged": False,
                           "tol": 1e-5 }
        self.fit(p0)
        self.params()

    def normal_form(self):
        raise NotImplementedError("The parameter fit does"
                                "not have an algebraic normal form")

    def sys(self, x, y, u, m):
        z = u[:2]
        a = u[2]
        b = u[3]
        alpha = u[4]
        phi = u[5:]

        # Convenience trig variables
        c = np.cos(phi)
        s = np.sin(phi)
        ca = np.cos(alpha)
        sa = np.sin(alpha)

        # function values (x0, y0, x1, y1, ....)
        f = flatzip(x.flatten() - z[0] - (a*ca*c - b*sa*s),
                    y.flatten() - z[1] - (a*sa*c + b*ca*s))
        # Jacobian
        J0 = np.zeros((2*m, m))
        for i in range(m):
            J0[2*i, i] = a*ca*s[i] + b*sa*c[i]
            J0[2*i+1, i] = a*sa*s[i] - b*ca*c[i]

        J = np.hstack((flatzip(-np.ones(m), np.zeros(m)).reshape(2*m, -1),
                       flatzip(np.zeros(m), -np.ones(m)).reshape(2*m, -1),
                       flatzip(-ca*c, -sa*c).reshape(2*m, -1),
                       flatzip(sa*s, -ca*s).reshape(2*m, -1),
                       flatzip(a*sa*c+ b*ca*s, -a*ca*c+b*sa*s).reshape(2*m, -1),
                       J0 ))
        return f, J

    def _check_radius(self, a, b):
        if abs(a-b)/(a+b) < self.properties["circ_tol"]:
            msg = "Ellipse is near-circular - nonlinear fit may not succeed"
            raise RuntimeError()

    def estimate_p0(self, x, y, p0):
        # p0 = (centerx, centery, a, b, phi)
        psi = np.arctan2(y-p0[1], x-p0[0]) - p0[4]
        return np.hstack([p0, psi.flatten()]).T

    def fit(self, p0):
        converged = False
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        u = self.estimate_p0(x, y, p0)

        # Iterate using Gauss Newton
        for nIts in range( self.properties['max_iterations'] ):
            f, J = self.sys(x, y, u, m)
            h = np.linalg.lstsq( -J, f ) [0]
            u = u + h
            # Check for convergence
            delta = np.linalg.norm(h, np.inf)/np.linalg.norm(u, np.inf)
            if delta < self.properties['tol']:
                converged = True
                break
        self.pfinal = u
        return tuple(u[:5]) + (converged, )

    def params(self):
        self.fdata.center = self.pfinal[:2]
        self.fdata.a = self.pfinal[2]
        self.fdata.b = self.pfinal[3]
        self.fdata.phi0 = self.pfinal[4]
        if self.fdata.a > self.fdata.b:
            eps = np.sqrt(1.0-(self.fdata.b/self.fdata.a)**2)
        else:
            eps = np.sqrt(1.0-(self.fdata.a/self.fdata.b)**2)
        self.fdata.eps = eps


class EllipseLevenberg(EllipseGaussNewton):

    def sys(self, u, x, y):
        # Convenience trig variables
        ca = np.cos(u[4])
        sa = np.sin(u[4])
        c = np.cos(u[5:])
        s = np.sin(u[5:])
        # function values (x0, y0, x1, y1, ....)
        f = flatzip(x.flatten() - u[0] - (u[2]*ca*c - u[3]*sa*s),
                    y.flatten() - u[1] - (u[2]*sa*c + u[3]*ca*s))
        return f

    def chi(self, x, y):
        return lambda p: self.sys(p, x, y)

    def jacobian(self, u, x, y):
        m = x.size
        c = np.cos(u[5:]) # cos(phi)
        s = np.sin(u[5:]) # sin(phi)
        ca = np.cos(u[4]) # cos(alpha)
        sa = np.sin(u[4]) # sin(alpha)

        J0 = np.zeros((2*m, m))
        for i in range(m):
            J0[2*i, i] = u[2]*ca*s[i] + u[3]*sa*c[i]
            J0[2*i+1, i] = u[2]*sa*s[i] - u[3]*ca*c[i]
        J = np.hstack((flatzip(-np.ones(m), np.zeros(m)).reshape(2*m, -1),
                       flatzip(np.zeros(m), -np.ones(m)).reshape(2*m, -1),
                       flatzip(-ca*c, -sa*c).reshape(2*m, -1),
                       flatzip(sa*s, -ca*s).reshape(2*m, -1),
                       flatzip(u[2]*sa*c+ u[3]*ca*s,
                               -u[2]*ca*c+u[3]*sa*s).reshape(2*m, -1),
                       J0 ))
        return J

    def fit(self, p0):
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        u = self.estimate_p0(x, y, p0)
        efunc = self.chi(x, y)
        jac = lambda p: self.jacobian(p, x, y)
        self.pfinal, err = leastsq(efunc, u, Dfun=jac)
        return tuple(u[:5]) + (err, )
