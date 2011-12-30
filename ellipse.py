# -*- coding: utf-8 -*-
"""
ellipse.py

Matplotlib Ellipse patch and allgebraic fit

"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'
__svn_id__ = '$Id$'

import numpy as np
import numpy.linalg as la
from collections import namedtuple
from matplotlib import patches
from matplotlib import lines
from matplotlib.pyplot import arrow
from scipy.optimize import leastsq

FIELDS = ['xraw', 'yraw', 'a', 'b',
          'eps', 'center', 'nform',
          'popt', 'paxes', 'phi0']

sq2 = np.sqrt(2)

def pol2cat(phi, rad, deg=True):
    if deg:
        phi = np.radians(phi)
    return rad * np.cos(phi), rad * np.sin(phi)

def ellipse_polar(phi, b, e, t):
    return  b/np.sqrt(1-(e*np.cos(phi-t))**2)

def rotmat(phi, inverse=False):
    if inverse:
        return np.matrix([[cos(phi),  sin(phi)],
                          [-sin(phi), cos(phi)]])
    else:
        return np.matrix([[cos(phi), -sin(phi)],
                          [sin(phi), cos(phi)]])


class EllipsePatch(patches.Ellipse):

    def __init__(self, xy, width, height, angle=0.0, **kwargs):
        if width < height:
            width, height = height, width
        super(EllipsePatch, self).__init__(xy, width, height, angle, **kwargs)
        self.angle = angle
        self.width = width
        self.height = height
        self.center = xy

    def paxes(self, *args, **kwargs):
        cx, cy = self.center
        sc = (self.width + self.height)*0.5
        ssin = np.sin(np.radians(self.angle))*sc
        scos = np.cos(np.radians(self.angle))*sc
        xline = lines.Line2D( (-scos+cx, scos+cx),
                              (-ssin+cy, ssin+cy), *args, **kwargs)
        yline = lines.Line2D( (-ssin+cx, ssin+cx),
                              (scos+cy, -scos+cy), *args, **kwargs)
        return xline, yline

    def text(self):
        return ('Parameters:\n'
                'cx: %.3f\n'
                'cy: %.3f\n'
                'a: %.3f\n'
                'b: %.3f\neps: %.3f'
                %(self.center[0],
                  self.center[1],
                  self.width/2.0,
                  self.height/2.0,
                  self.eps))

    @property
    def eps(self):
        if self.width < self.height:
            eps = np.sqrt(1- self.width**2/self.height**2)
        else:
            eps = np.sqrt(1- self.height**2/self.width**2)
        return eps

class EllipseBase(object):

    def __init__(self, x, y):
        self.fdata = namedtuple('fdata', FIELDS )
        self.fdata.xraw = x
        self.fdata.yraw = y

    def __getattr__(self, attr):
        if hasattr(self.fdata, attr) and not attr.startswith('_'):
            return getattr(self.fdata, attr)
        else:
            return getattr(self, attr)

    def normal_form(self):
        """Transforms an arbitrary conic equation into its normalform"""
        par = self.fdata.popt
        A = np.matrix([[par[0], par[1]/sq2], [par[1]/sq2, par[2]]])
        val, vec = la.eig(A)
        b = np.matrix(par[3:5]).T
        t = (-1./2.*b.T*A.I).T
        c = t.T*A*t + b.T*t + par[-1]
        self.fdata.center = np.array(t).flatten()
        self.fdata.nform = np.append(val, c)
        self.fdata.paxes = vec
        if c > 0.:
            self.fdata.nform *= -1.0
        return self.fdata.nform

    def params(self):
        """Transforms the normal form into its parameter representation"""
        nform = self.fdata.nform
        pa = self.fdata.paxes
        a = np.sqrt(-nform[2]/nform[0])
        b = np.sqrt(-nform[2]/nform[1])
        if a > b:
            eps = np.sqrt(1.0-(b/a)**2)
            phi0 = np.arctan(pa[1, 0]/pa[0, 0])
        else:
            eps = np.sqrt(1.0-(a/b)**2)
            phi0 = np.arctan(pa[1, 1]/pa[0, 1])
        self.fdata.a = a
        self.fdata.b = b
        self.fdata.eps = eps
        self.fdata.phi0 = phi0

class Ellipse(EllipseBase):
    """Fit cartesian data to an ellipse without any constraints"""

    def __init__(self, x, y):
        super(Ellipse, self).__init__(x, y)
        self.fit()
        self.normal_form()
        self.params()

    def fit(self):
        """Returns the design matrix for an algebraic fit of an ellipsis"""
        m = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        SC = np.matrix(np.hstack([x**2, sq2*x*y, y**2,
                                  x, y, np.ones((m, 1))]))
        # design matrix
        DM = SC.T*SC
        w, v = la.eig(DM)
        w, v = self.lambda_min(w, v)
        self.fdata.popt = v.flatten()

    def lambda_min(self, w, v):
        """Returns the minimum eigenvalue and the corresponding eigenvector"""
        ind = np.where(w==w.min())[0]
        return np.array(w[ind]), np.array( v[:, ind])

class EllipseBookstein(EllipseBase):

    def __init__(self, x, y):
        super(EllipseBookstein, self).__init__(x, y)
        self.fit()
        self.normal_form()
        self.params()

    def fit(self):
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)

        # B *[v; w] = 0, with the constraint norm(w) == 1
        B = np.hstack([ x, y, np.ones((m, 1)), x**2, sq2*x*y, y**2 ])

        # To enforce the constraint, we need to take the QR decomposition
        Q, R = la.qr(B)
        # Decompose R into blocks
        R11 = R[0:3, 0:3]
        R12 = R[0:3, 3:6]
        R22 = R[3:6, 3:6]

        # Solve R22 * w = 0 subject to norm(w) == 1
        U, S, V = la.svd(R22)
        V = V.T
        w = V[:, 2]
        # Solve for the remaining variables
        v = np.dot(la.solve( -R11, R12 ), w)
        self.fdata.popt = np.append(w, v).flatten()

class EllipseTrace(EllipseBase):

    def __init__(self, x, y):
        super(EllipseTrace, self).__init__(x, y)
        self.fit()
        self.normal_form()
        self.params()

    def fit(self):
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)

        # Coefficient matrix
        B = np.hstack([ sq2* x* y, y**2-x**2, x, y, np.ones((m, 1)) ])
        B = np.matrix(B)
        v = la.lstsq(B, -x**2)[0]
        self.fdata.popt = np.append([1-v[1], v[0], v[1]], v[2:5])
