"""
ellipse.py

Matplotlib EllipsePatch and algebraic and geometric
fit routines.

"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'
__svn_id__ = '$Id$'

import numpy as np
import numpy.linalg as la
from scipy.optimize import leastsq
from collections import namedtuple
from matplotlib import patches
from matplotlib import lines

FIELDS = ['xraw', 'yraw', 'a', 'b',
          'eps', 'center', 'nform',
          'popt', 'paxes', 'phi0']

sq2 = np.sqrt(2)

def pol2cat(phi, rad, deg=True):
    """Converts polar coordinates to cartesian coordinates
    using numpy ufuncs"""
    if deg:
        phi = np.radians(phi)
    return rad * np.cos(phi), rad * np.sin(phi)

def ellipse_polar(phi, b, eps, alpha):
    """Return r(phi) of an ellipse providing the parameters:
    b (semiaxis), eps (excentricity), alpha (rotation)
    """
    return  b/np.sqrt(1-(eps*np.cos(phi-alpha))**2)

def flatzip(x, y):
    """Zipping two ndarrays and flatten the result"""
    zipped = np.array(zip(x, y))
    return zipped.flatten()


class EllipsePatch(patches.Ellipse):

    def __init__(self, center, width, height, angle=0.0, **kwargs):
        if width < height:
            width, height = height, width
        super(EllipsePatch, self).__init__(center, width, height,
                                           angle, **kwargs)
        self.angle = angle
        self.width = width
        self.height = height
        self.center = center

    def paxes(self, *args, **kwargs):
        """Returns the primary axes as matplotlib.lines.Line2D patches"""
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
        """Fitted parameters as string"""
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
        """Returns the excentricity of the ellipse"""
        if self.width < self.height:
            eps = np.sqrt(1- self.width**2/self.height**2)
        else:
            eps = np.sqrt(1- self.height**2/self.width**2)
        return eps

class EllipseBase(object):

    def __init__(self, x, y):
        """x and y is the raw data to fit"""
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

    def normal_form(self):
        """Transforms an arbitrary conic equation into its normal form"""
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
        """Transforms the normal form into its parameter representation."""
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
        """Fits the data"""
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
    """Fit an ellipse with the Bookstein constraint

    lamdba1**2 + lambda2**2 = 1
    where lambda are the eigenvalues of the matrix A
    """

    def __init__(self, x, y):
        super(EllipseBookstein, self).__init__(x, y)
        self.fit()
        self.normal_form()
        self.params()

    def fit(self):
        """Fits the raw data"""
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
    """Fit an ellipse with the trace constraint

    lamdba1 + lambda2 = 1
    where lambda are the eigenvalues of the matrix A
    """

    def __init__(self, x, y):
        super(EllipseTrace, self).__init__(x, y)
        self.fit()
        self.normal_form()
        self.params()

    def fit(self):
        """Fit the data"""
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)

        # Coefficient matrix
        B = np.hstack([ sq2* x* y, y**2-x**2, x, y, np.ones((m, 1)) ])
        B = np.matrix(B)
        v = la.lstsq(B, -x**2)[0]
        self.fdata.popt = np.append([1-v[1], v[0], v[1]], v[2:5])

class EllipseGaussNewton(EllipseBase):
    """Geometric fit (best fit) using a Gauss-Newton Iteration"""

    def __init__(self, x, y, p0):
        """x and y is the raw data, p0 is the initial estimate of the
        ellipse parameters (in polar coordinates)
        """
        super(EllipseGaussNewton, self).__init__(x, y)
        self.properties = {"circ_tol": 1e-5,
                           "max_iterations": 666,
                           "converged": False,
                           "tol": 1e-5 }
        self.fit(p0)
        self.params()

    def normal_form(self):
        """Raise a NotImplementedError since the normal form
        does not calculated. It just overwrites the method inherited
        """
        raise NotImplementedError("The parameter fit does"
                                "not have an algebraic normal form")

    def parfunc(self, x, y, p):
        """Returns f(x, y, p), J(x, y, p), where f are the function values
        and J is the Jacobian."""
        m = x.size
        z = p[:2]
        a = p[2]
        b = p[3]
        alpha = p[4]
        phi = p[5:]

        # Convenience trigometric variables
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
        """Consider ellipse as circle if excentricity is
        smaller than properties['circ_tol']"""
        if abs(a-b)/(a+b) < self.properties["circ_tol"]:
            msg = "Ellipse is near-circular - nonlinear fit may not succeed"
            raise RuntimeError()

    def estimate_p0(self, x, y, p0):
        """Esitmates start paramters for non-linear leastsq

        Meaning of p0:
        p0 = (centerx, centery, a, b, phi)
        """
        psi = np.arctan2(y-p0[1], x-p0[0]) - p0[4]
        return np.hstack([p0, psi.flatten()]).T

    def fit(self, p0):
        """Fit the data"""
        converged = False
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        p = self.estimate_p0(x, y, p0)

        # Iterate using Gauss Newton
        for nIts in range( self.properties['max_iterations'] ):
            f, J = self.parfunc(x, y, p)
            h = np.linalg.lstsq( -J, f ) [0]
            p = p + h
            # Check for convergence
            delta = np.linalg.norm(h, np.inf)/np.linalg.norm(p, np.inf)
            if delta < self.properties['tol']:
                converged = True
                break
        self.pfinal = p
        return tuple(p[:5]) + (converged, )

    def params(self):
        """Slice out ellipse parameters out of the final fit parameters"""
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
    """Geometric fit ("best fit") using the Levenberg-Marquard algorithm"""

    def parfunc(self, p, x, y):
        """Returns f(x, y, p), the function values."""

        # Convenience trig variables
        ca = np.cos(p[4])
        sa = np.sin(p[4])
        c = np.cos(p[5:])
        s = np.sin(p[5:])
        # function values (x0, y0, x1, y1, ....)
        f = flatzip(x.flatten() - p[0] - (p[2]*ca*c - p[3]*sa*s),
                    y.flatten() - p[1] - (p[2]*sa*c + p[3]*ca*s))
        return f

    def chi(self, x, y):
        return lambda p: self.parfunc(p, x, y)

    def jacobian(self, p, x, y):
        """Returns jacobian"""
        m = x.size
        c = np.cos(p[5:]) # cos(phi)
        s = np.sin(p[5:]) # sin(phi)
        ca = np.cos(p[4]) # cos(alpha)
        sa = np.sin(p[4]) # sin(alpha)

        J0 = np.zeros((2*m, m))
        for i in range(m):
            J0[2*i, i] = p[2]*ca*s[i] + p[3]*sa*c[i]
            J0[2*i+1, i] = p[2]*sa*s[i] - p[3]*ca*c[i]
        J = np.hstack((flatzip(-np.ones(m), np.zeros(m)).reshape(2*m, -1),
                       flatzip(np.zeros(m), -np.ones(m)).reshape(2*m, -1),
                       flatzip(-ca*c, -sa*c).reshape(2*m, -1),
                       flatzip(sa*s, -ca*s).reshape(2*m, -1),
                       flatzip(p[2]*sa*c+ p[3]*ca*s,
                               -p[2]*ca*c+p[3]*sa*s).reshape(2*m, -1),
                       J0 ))
        return J

    def fit(self, p0):
        """Fit the data"""
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        u = self.estimate_p0(x, y, p0)
        efunc = self.chi(x, y)
        jac = lambda p: self.jacobian(p, x, y)
        self.pfinal, err = leastsq(efunc, u, Dfun=jac)
        return tuple(self.pfinal[:5]) + (err, )
