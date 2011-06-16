#!/usr/bin/python

from numpy import ones, cos, matrix, pi, arange, degrees, zeros, append
from numpy import degrees, radians, array, round, cos, sin, sqrt
from numpy.random import randn, rand
from numpy import where
from pylab import plot, show, xlim, figure,title, polar, gca
from pdb import set_trace
from numpy import linalg as LA
from scipy.optimize import curve_fit


sens = 1.0
nmeas = 100.0
eamp = 0.01
err = eamp*randn(nmeas)
r0 = ones(nmeas)*10
phi = arange(nmeas)/nmeas*2*pi  # angles of the channels


def ellipse_polar(phi, b, e, t):    

    rad = b/sqrt(1-(e*cos(phi-t))**2)
    
    return rad


def pol2cat(rad, phi, deg=True):

    if deg:
        phi = radians(phi)

    return rad * cos(phi), rad * sin(phi)


def scatter_matrix(points):
    """Returns the design matrix for an algebraic fit of an ellipsis
    """
    
    DM = array([])
    for (x, y) in points:
        DM = append(DM, array([x**2, 2*x*y, y**2, x, y, 1]))
    DM= matrix(DM.reshape(-1, 6))

    return DM.T*DM

def est_pinit(xdata, ydata):
    """Estimates inital fit paramters"""
    
    
    ind = where(xdata==xdata.max())
    phi_max = phi[ind]
    
    return array([ydata.mean(), 0.5, phi_max])


pfunc = array([5, 0.5, rand()*pi])
r0 = ellipse_polar(phi, *pfunc)

# non linear lsq
pinit = est_pinit(r0+err, phi)
popt, pcov = curve_fit(ellipse_polar, phi, r0+err, p0=pinit)
print popt
print pfunc

x, y = pol2cat(r0+err, phi)
points = [(i, j) for i, j in zip(x, y)]

S = scatter_matrix(points)
w, v = LA.eig(S)



# x = matrix(err)*A.I
figure(1)
#polar(phi, r0, '-r', lw=1.5)
polar(phi, ellipse_polar(phi, *popt), '-r', lw=1.5)
#polar(phi, r0, '-b', lw=1.5)
polar(phi, r0 + err, 'ob')
#gca().set_rmax(r0.max()*1.05)

# set_trace()

#polar(phi, d0 + x, 'or')
title('Deviation')

show()



