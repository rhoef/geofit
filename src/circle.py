"""
circle.py

Matplotlib CirclePatch and algebraic and geometric
fit routines.
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'

import numpy as np
import numpy.linalg as la
from scipy.optimize import leastsq
from matplotlib import patches
from collections import namedtuple

FIELDS = ['xraw', 'yraw', 'center', 'popt', 'radius']

class CirclePatch(patches.Circle):

    def __init__(self, *args, **kw):
        super(CirclePatch, self).__init__(*args, **kw)

class CircleBase(object):

    def __init__(self, x, y):
        """x and y is the raw data to fit."""
        self.fdata = namedtuple('fdata', FIELDS )
        self.fdata.xraw = x
        self.fdata.yraw = y

    def fit(self):
        """This method needs to be overwritten in child classes"""
        raise NotImplementedError

    def __getattr__(self, attr):
        """To treat attributes from fdata as attributes from self"""
        if hasattr(self.fdata, attr) and not attr.startswith('_'):
            return getattr(self.fdata, attr)
        else:
            return getattr(self, attr)

    def text(self):
        """Fitted parameters as string"""
        return ('Parameters:\n'
                'cx: %.3f\n'
                'cy: %.3f\n'
                'r: %.3f\n'
                %(self.center[0], self.center[1], self.radius))

class CircleAlgebraic(CircleBase):
    """Algebraic fit cartesian data to a circle

    Access the data and fit parameters:
    >>>circle = CircleAlgebraic(x, y)
    >>>c = circle.center
    >>>r = circle.radius

    Constants of the algebraic equation are stored in
    >>>circle.popts
    """

    def __init__(self, x, y):
        super(CircleAlgebraic, self).__init__(x, y)
        self.fit()
        self.params()

    def fit(self):
        """Fits the data"""
        m = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        B = np.matrix(np.hstack([x**2 + y**2, x, y, np.ones((m, 1))]))

        # design matrix
        DM = B.T*B
        w, v = la.eig(DM)
        w, v = self.lambda_min(w, v)
        self.fdata.popt = v.flatten()

    def lambda_min(self, w, v):
        """Returns the minimum eigenvalue and the corresponding eigenvector"""
        ind = np.argmin(w)
        return np.array(w[ind]), np.array( v[:, ind])

    def params(self):
        """Transforms the normal form into its parameter representation."""
        a = self.fdata.popt[0]
        b = self.fdata.popt[1:3]
        c = self.fdata.popt[3]

        self.fdata.radius = np.sqrt(np.linalg.norm(b, 2)**2/(4*a**2) - c/a)
        self.fdata.center = np.array([-b[0]/(2*a), -b[1]/(2*a)])

class CircleLevenberg(CircleAlgebraic):
    """Fit cartesian data to an cirlce using the Levenberg-Marquardt"""

    def __init__(self, x, y, p0):
        super(CircleAlgebraic, self).__init__(x, y)
        self.p0 = p0
        self.fit()
        self.params()

    def lambda_min(self):
        raise NotImplementedError

    def func(self, p, x, y):
        return  (np.sqrt((x-p[0])**2 + (y-p[1])**2) - p[2]).flatten()

    def chi(self, x, y):
        return lambda p: self.func(p, x, y)

    def jacobian(self, p, x, y):
        """Jacobian of the circle equation"""
        m = x.size
        denom = np.sqrt((p[0]-x)**2 + (p[1]-y)**2)
        J = np.hstack([(p[0]-x)/denom, (p[1]-y)/denom, -np.ones((m, 1))])
        return J

    def params(self):
        self.fdata.center = self.fdata.popt[:2]
        self.fdata.radius = self.fdata.popt[2]

    def fit(self):
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        efunc = self.chi(x, y)
        jac = lambda p: self.jacobian(p, x, y)
        pfinal, err = leastsq(efunc, self.p0, Dfun=jac)
        self.fdata.popt = pfinal
        return tuple(pfinal) + (err, )
