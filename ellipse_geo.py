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
from ellipse import rotmat
import numpy as np
from scipy.optimize import leastsq


class EllipseGeometric(EllipseBase):

    def __init__(self, x, y, p0):
        super(EllipseGeometric, self).__init__(x, y)
        self.params = {"circ_tol": 1e-5,
                       "max_iterations": 666,
                       "converged": False}
        self.fit(p0)
#        self.params()

    def normal_form(self):
        raise NotImplementError("The parameter fit does"
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

        # Rotation matrices
        Q    = np.array( [[ca, -sa], [sa, ca]] )
        Qdot = np.array( [[-sa, -ca], [ca, -sa]] )


        f = np.vstack((x[0, :] - z[0] - (a*ca*c - b*sa*s),
                       x[1, :] - z[1] - (a*sa*c + b*ca*s))).T.flatten()

        ones = np.ones((m, 1))
        zeros = np.zeros((m, 1))
        J = np.hstack((-ones, zeros,
                        zeros, ones,
                        -(ca*c).reshape((m, 1)), -(sa*c).reshape((m, 1)),
                        (sa*s).reshape((m, 1)), -(ca*s).reshape((m, 1)),
                        (a*sa*c+ b*ca*s).reshape((m, 1)),
                        (-a*ca*c+b*sa*s).reshape((m, 1)),
                        np.zeros((m, 2*m))
                        ))
        for i in range(m):
            J[i, 2*i+5] = a*ca*s[i] + b*sa*c[i]
            J[i, 2*i+6] = a*sa*s[i] - b*ca*c[i]
        import pdb; pdb.set_trace()
        return f, J

    def _check_radius(self, a, b):
        if abs(a-b)/(a+b) < self.params["circ_tol"]:
            msg = "Ellipse is near-circular - nonlinear fit may not succeed"
            raise RuntimeError()

    def estimate_p0(self, x, y, p0):
        # p0 = (centerx, centery, a, b, phi)
        psi = np.arctan2(y-p0[1], x-p0[0]) - p0[4]
        return np.hstack([p0, psi.flatten()]).T

    def fit(self, p0):
        m  = self.fdata.xraw.size
        x = self.fdata.xraw.reshape(m, 1)
        y = self.fdata.yraw.reshape(m, 1)
        u = self.estimate_p0(x, y, p0)

        # Iterate using Gauss Newton
        for nIts in range( self.params['max_iterations'] ):
            # Find the function and Jacobian
            f, J = self.sys(x, y, u, m)
            # Solve for the step and update u
            # h = linalg.solve( -J, f )
            h = np.linalg.lstsq( -J, f ) [0]
            u = u + h

            # Check for convergence
            delta = np.linalg.norm(h, inf)/np.linalg.norm(u, inf)
            if delta < params['tol']:
                converged = True
                break

        # alpha = u[-5]
        # a     = u[-4]
        # b     = u[-3]
        # z     = u[-2:]
        self.pfinal = u
        # return  (centerX, centerY, a, b, phi)
        return tuple(u[:5]) + (converged, )
