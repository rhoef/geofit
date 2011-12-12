#!/usr/bin/python

from matplotlib import use
use("TkAgg")

import numpy as np
import pylab as pl
from numpy.random import randn, rand
from matplotlib.pyplot import figure, show, draw, gca

from pdb import set_trace
from numpy import linalg as LA
from matplotlib.patches import Circle

from ellipse import Ellipse
from ellipse import EllipseGeometric
from ellipse import EllipsePatch
from ellipse import ellipse_polar
from ellipse import pol2cat

def fit_ellipse(x, y):

    # fit raw data
    el  = Ellipse(x, y)
    elg = EllipseGeometric(x, y, p0=(el.center[0], el.center[1],
                                     el.a, el.b, el.phi0))
    # plot fit in normal form
    el_hl1 = EllipsePatch( (0, 0), 2*el.a, 2*el.b, fill=False, color='g')
    # plot fit of raw data
    el_raw = EllipsePatch( tuple(el.center), 2*el.a, 2*el.b,
                           angle=np.degrees(el.phi0), fill=False, color='r',
                           lw=2)

    el_gfit = EllipsePatch( tuple(elg.center), 2*elg.a, 2*el.b,
                            angle=np.degrees(el.phi0), fill=False, color='y',
                            lw=2)


    fig, ax = getAxes()

    ax.plot(x, y, 'bo')
    ax.add_artist(el_hl1)
    ax.add_artist(el_gfit)
    ax.add_artist(el_raw)
    xline, yline = el_raw.paxes(ls='--', color='k')
    ax.add_artist(xline)
    ax.add_artist(yline)
    fig.text(0.05, 0.40, el_raw.text(), fontsize=11)
    setLimits(ax)

    print 'Center: %.3f, %.3f' %tuple(el.center)
    print 'a, b: %.3f, %.3f' %(el.a, el.b)
    print 'Alpha: %.3f' %el.phi0
    print 'Center: %.3f, %.3f' %tuple(elg.center)
    print 'a, b: %.3f, %.3f' %(elg.a, elg.b)
    print 'Alpha: %.3f' %elg.phi0


def setLimits(ax):
    lmax = max(list(ax.get_xlim())+list(ax.get_ylim()))
    ax.set_ylim((-lmax, lmax))
    ax.set_xlim(ax.get_ylim())
    ax.set_title('Deviation')

def getAxes():
    fig = figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.axhline(y=0.0, color='k')
    ax.axvline(x=0.0, color='k')
    return fig, ax


if __name__ == '__main__':

    sens = 1.0
    nmeas = 100
    err = randn(nmeas)*0.1
    phi = np.linspace(0, 2*np.pi, nmeas)  # angles of the channels
    r0 = ellipse_polar(phi, 4, rand(), rand()*np.pi)
    x, y = pol2cat(phi, r0+err, deg=False)

    # data= np.loadtxt('dipoletest2.csv', delimiter=',')
    # x = data[1:, 0] - data[0, 0]
    # y = data[1:, 1] - data[0, 1]
#    import pdb; pdb.set_trace()

    el = fit_ellipse(x, y)



    show()

