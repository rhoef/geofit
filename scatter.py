#!/usr/bin/python

from numpy import ones, cos, matrix, pi, arange, degrees, zeros, append, linspace
from numpy import degrees, radians, array, round, cos, sin, sqrt, cov
from numpy.random import randn, rand
from numpy import where
from pylab import plot, show, xlim, figure,title, polar, gca, ylim
from pdb import set_trace
from numpy import linalg as LA
from matplotlib.path import Path


sens = 1.0
nmeas = 5000.0
eamp = 0.5
err = eamp*randn(nmeas)
phi = arange(nmeas)/nmeas*2*pi  # angles of the channels



def ellipse_polar(phi, b, e, t):    

    rad = b/sqrt(1-(e*cos(phi-t))**2)
    
    return rad


def xlimits(a, b, c, d, e, f):
    
    p = (b*e - 2*d*c)/(2*(b**2 - c * a))
    q = (e**2 - 4 * c * f)/(4*(b**2 - c * a))
    
    xmin = -p/2 - sqrt((p/2)**2 - q)
    xmax = -p/2 + sqrt((p/2)**2 - q)

    return xmin*0.999, xmax*0.999


def ellipse(par, n=300):

    (a, b, c, d, e, f) = par.flatten()
    xlim = xlimits(a, b, c, d, e, f)
    x = linspace(xlim[0], xlim[1], n)
    
    p = (2*b*x+e)/c
    q = (a*x**2 + d*x +f)/c
    
    y = -p/2 + sqrt((p/2)**2 - q) # upper half
    y = append(y, (-p/2 - sqrt((p/2)**2 - q))[::-1]) # lower half

    x = append(x, x[::-1])

    return x, y

def lambda_min(w, v):
    
    ind = where(w==w.min())[0]
    
    return array(w[ind]), array( v[:, ind])   


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

    return DM

r0 = ellipse_polar(phi, 10, rand(), rand()*pi)

x, y = pol2cat(r0+err, phi, deg=False)

#COV = cov(x, y)
#pl, pa = LA.eig(COV)

#set_trace()


points = [(i, j) for i, j in zip(x, y)]
S = scatter_matrix(points)
DM = S.T*S
w, v = LA.eig(DM)
w, v = lambda_min(w, v)

xf, yf = ellipse(v)

#set_trace()
#COV = cov(xf, yf)
#pl, pa = LA.eig(COV)

#COV2 = cov(x, y)
#pl2, pa2 = LA.eig(COV2)

set_trace()

# x = matrix(err)*A.I
figure(1)
plot(x, y, 'bo')
plot(xf, yf, '-r')

ylim(xlim())
sc =100
#plot(array([0, pa[0,1]])*sc, array([0, pa[0,0]])*sc, '-y', lw=2)
#plot(array([0, pa[1,1]])*sc, array([0, pa[1,0]])*sc, '-g', lw=2)

#plot(array([0, pa[0,0]])*sc, array([0, pa[1,0]])*sc, '-y', lw=2)
#plot(array([0, pa[0,1]])*sc, array([0, pa[1,1]])*sc, '-g', lw=2)
#plot(array([0, pa2[0,0]])*sc, array([0, pa2[1,0]])*sc, '-r', lw=2)
#plot(array([0, pa2[0,1]])*sc, array([0, pa2[1,1]])*sc, '-k', lw=2)


# set_trace()
#print pa
#print pl
# polar(phi, ellipse_polar(phi, *popt), '-r', lw=1.5)
#polar(phi, r0, '-b', lw=1.5)
# polar(phi, r0 + err, 'ob')
# gca().set_rmax(r0.max()*1.05)

# set_trace()

# polar(phi, d0 + x, 'or')
title('Deviation')

show()



